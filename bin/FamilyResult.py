"""Sequence/AGP/provenance publication for a completed family search."""
from dataclasses import asdict
from html import escape
import os
from pathlib import Path
import tempfile

from ExtensionFamilies import atomic_json, fingerprint, read_json
from FamilyOverlap import joined_sequence, read_sequence
from FamilyTools import content_hash

SCHEMA = "degap-family-result-v1"


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def load_result(directory):
    """Load only a complete, content-verified publication, never log heuristics."""
    directory = Path(directory)
    result = read_json(directory / "family-result.json")
    payload = {k: v for k, v in result.items() if k != "payload_sha256"}
    if result.get("schema") != SCHEMA or result.get("payload_sha256") != fingerprint(payload):
        raise ValueError("Invalid extension-family result manifest")
    if (result.get("status") == "closed") != (result.get("winner") is not None):
        raise ValueError("Family result closure state is inconsistent")
    for relative, digest in result["artifacts"].items():
        path = directory / relative
        if path.resolve().parent != directory.resolve() or not path.is_file() or content_hash(path) != digest:
            raise ValueError(f"Family result artifact changed: {path}")
    for required in (result["final_sequence"], "family-components.fa", "process.agp", "process.log", "process.summary", "family-report.html"):
        if required not in result["artifacts"]:
            raise ValueError(f"Family result is missing required artifact: {required}")
    return result


def publish_search(search, name=None):
    """Export one joined sequence, or separate unclosed family sequences.

    AGP components are actual exported cumulative sequences, not invented reads
    or falsely attributed original contigs. Family/source provenance accompanies
    them in the manifest. Standard AGP coordinates are one-based inclusive.
    """
    if search.status == "running" or search.pending:
        raise ValueError("Cannot publish an incomplete family search")
    directory = search.directory
    name = name or directory.name
    if not name or Path(name).name != name or any(c.isspace() for c in name):
        raise ValueError("Result name must be a single whitespace-free filename component")
    result_file = directory / "family-result.json"
    if result_file.exists():
        result = load_result(directory)
        if result["search_contract"] != search.contract_id or result["status"] != search.status:
            raise ValueError("Published result belongs to a different search")
        return result

    final_file = name + ".final.fa"
    result = {
        "schema": SCHEMA, "search_contract": search.contract_id,
        "status": search.status, "total_rounds": search.total_rounds,
        "total_extension_length": search.total_length,
        "final_sequence": final_file, "winner": asdict(search.winner) if search.winner else None,
        "gap_length": None, "families": [], "components": {}, "artifacts": {},
    }
    for family in search.families.values():
        result["families"].append({
            "key": family.key, "direction": family.direction, "attempts": family.attempts,
            "added_bases": family.added_bases, "stop_reason": family.stop_reason,
            "versions": [asdict(v) for v in family.versions],
        })
    components = []
    agp = []

    def component(version):
        version.verify()
        key = f"family{search.order[version.family]}.v{version.number}"
        if key not in result["components"]:
            components.append(f">{key}\n{read_sequence(version.path)}\n")
            result["components"][key] = asdict(version)
        return key

    if search.winner:
        winner = search.winner
        left_family, right_family = search.families[winner.left_family], search.families[winner.right_family]
        left, right = left_family.versions[winner.left_version], right_family.versions[winner.right_version]
        sequence = joined_sequence(winner, left, right)
        final_text = f">{name}\n{sequence}\n"
        result["gap_length"] = len(sequence) - left_family.versions[0].length - right_family.versions[0].length
        left_id, right_id = component(left), component(right)
        agp.append(f"{name}\t1\t{winner.left_end}\t1\tW\t{left_id}\t1\t{winner.left_end}\t+")
        if right.length > winner.right_end:
            agp.append(f"{name}\t{winner.left_end + 1}\t{len(sequence)}\t2\tW\t{right_id}\t{winner.right_end + 1}\t{right.length}\t+")
    else:
        records = []
        for family in search.families.values():
            version = family.current
            key = component(version)
            record_id = f"{name}.family{search.order[family.key]}.{search.status}"
            records.append(f">{record_id}\n{read_sequence(version.path)}\n")
            agp.append(f"{record_id}\t1\t{version.length}\t1\tW\t{key}\t1\t{version.length}\t+")
        final_text = "".join(records)

    summary = ["family\tdirection\tattempts\tversions\tadded_bases\tstop_reason"]
    for family in search.families.values():
        summary.append("\t".join(map(str, (family.key, family.direction, family.attempts,
                                           len(family.versions), family.added_bases, family.stop_reason or ""))))
    log = ["DEGAP extension-family search", f"Status: {search.status}",
           f"Total attempted rounds: {search.total_rounds}", f"totalExtensionLength: {search.total_length}"]
    if search.winner:
        log += ["GAP can be closed!", f"GAP Length: {result['gap_length']}",
                f"Selected families: {search.winner.left_family}/{search.winner.left_version} + "
                f"{search.winner.right_family}/{search.winner.right_version}"]
    else:
        log += ["Not a closed gap: exported family sequences remain separate."]
    log.append(f"Final EXtendFile: {directory / final_file}")
    outputs = {
        final_file: final_text, "family-components.fa": "".join(components),
        "process.agp": "\n".join(agp) + "\n", "process.log": "\n".join(log) + "\n",
        "process.summary": "\n".join(summary) + "\n",
    }
    outputs["family-report.html"] = render_report(search, result)
    for relative, text in outputs.items():
        atomic_text(directory / relative, text)
        result["artifacts"][relative] = content_hash(directory / relative)
    result["payload_sha256"] = fingerprint(result)
    atomic_json(result_file, result)
    return result


def render_report(search, result):
    """Portable evidence tables for both chains; no single-direction attribution."""
    rows = []
    histories = []
    for family in search.families.values():
        values = (family.key, family.direction, family.attempts, family.added_bases, family.stop_reason or "—")
        rows.append("<tr>" + "".join(f"<td>{escape(str(v))}</td>" for v in values) + "</tr>")
        for version in family.versions:
            histories.append("<tr>" + "".join(f"<td>{escape(str(v))}</td>" for v in
                             (family.key, version.number, version.length, version.sha256)) + "</tr>")
    comparisons = []
    for candidates in search.checked.values():
        for candidate in candidates:
            selected = candidate == result["winner"]
            values = (f"{candidate['left_family']}/{candidate['left_version']}",
                      f"{candidate['right_family']}/{candidate['right_version']}",
                      candidate["identity"], candidate["overlap"], candidate["overhang"],
                      "Selected / 已采用" if selected else "Not selected / 未采用")
            comparisons.append("<tr>" + "".join(f"<td>{escape(str(v))}</td>" for v in values) + "</tr>")
    return f'''<!doctype html><html lang="en"><meta charset="utf-8">
<title>DEGAP family extension / 双向家族延伸</title>
<style>body{{font:16px system-ui;max-width:1100px;margin:32px auto;padding:0 16px;color:#243043}}
table{{border-collapse:collapse;width:100%;margin:16px 0}}th,td{{border:1px solid #cbd5e1;padding:8px;text-align:left;overflow-wrap:anywhere}}
th{{background:#edf2f7}}.scroll{{overflow:auto}}a{{color:#0369a1}}</style>
<h1>DEGAP family extension / 双向家族延伸</h1>
<p>Status / 状态: <strong>{escape(search.status)}</strong> · Attempts / 总尝试轮数: {search.total_rounds}
· Added bases / 总延伸长度: {search.total_length}</p>
<p>Only a closed result joins both families. Unclosed sequences are exported separately.<br>
只有闭合结果才连接两侧；未闭合时各家族序列分别导出。</p>
<p><a href="{escape(result['final_sequence'], quote=True)}">FASTA</a> · <a href="process.agp">AGP</a> ·
<a href="family-components.fa">AGP components / 组成序列</a> · <a href="family-result.json">Result provenance / 结果溯源</a> ·
<a href="families.json">History and checkpoint / 历史与断点</a></p>
<h2>Families / 延伸家族</h2><table><tr><th>Family</th><th>Direction</th><th>Attempts</th><th>Added bases</th><th>Stop reason</th></tr>{''.join(rows)}</table>
<h2>Eligible connections / 合格连接</h2><p>Rank: identity, overlap length, endpoint overhang, fixed version order.<br>
依次按一致性、重叠长度、端点未对齐长度和固定版本顺序选优；只采用一个结果。</p>
<table><tr><th>Left version</th><th>Right version</th><th>Identity %</th><th>Overlap</th><th>Overhang</th><th>Selection</th></tr>{''.join(comparisons)}</table>
<h2>Sequence history / 序列历史</h2><div class="scroll"><table><tr><th>Family</th><th>Version</th><th>Length</th><th>SHA-256</th></tr>{''.join(histories)}</table></div>
<p>Read mappings, local assemblies and their receipts remain under attempts/ for inspection.<br>
reads 比对、局部组装及其记录保存在 attempts/ 中，供复核。</p></html>'''
