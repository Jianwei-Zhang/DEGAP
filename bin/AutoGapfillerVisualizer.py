#!/usr/bin/env python3
"""
Visualizer - Step 5

生成 gap filling 结果的交互式 HTML 可视化报告。

功能：
1. 总览统计信息（总 gaps、填充率、BP 变化等）
2. 染色体视图（显示所有 gaps 的位置和填充状态）
3. Gap 详细信息（点击查看每个 gap 的详细报告）
4. 统计图表（填充类型分布、长度分布等）
5. 国际化支持（中英文切换）

输入：
- Gap filling 统计结果（gap_filling_results.json）
- Genome integration 结果（integration_config.json）
- Reads 统计信息（preprocessing_config.json）

输出：
- gap_filling_report.html：交互式可视化报告
"""

import json
import sys
from pathlib import Path
from datetime import datetime

try:
    from ReportStyle import esc, report_css
except ImportError:
    from .ReportStyle import esc, report_css


class GapFillingVisualizer:
    """生成 gap filling 结果的可视化报告"""

    def __init__(self, output_dir):
        """
        初始化 Visualizer

        Args:
            output_dir: 主输出目录
        """
        self.main_output_dir = Path(output_dir)

        # 输出目录
        self.output_dir = self.main_output_dir / "05.visualization"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 读取配置文件
        self._load_configurations()
    
    def _load_configurations(self):
        """读取所有需要的配置文件"""
        # Gap filling 结果
        gap_filling_results_file = self.main_output_dir / "03.gap_filling" / "gap_filling_results.json"
        if not gap_filling_results_file.exists():
            raise FileNotFoundError(
                f"Gap filling results not found: {gap_filling_results_file}\n"
                f"Please run AutoGapfillerPart3.py first!"
            )

        with open(gap_filling_results_file, 'r') as f:
            self.gap_filling_results = json.load(f)

        # Genome integration 结果（可选）
        integration_config_file = self.main_output_dir / "04.genome_integration" / "integration_config.json"
        if integration_config_file.exists():
            with open(integration_config_file, 'r') as f:
                self.integration_config = json.load(f)
        else:
            self.integration_config = None
            print("  Note: Genome integration results not found, some statistics may be unavailable")

        # Reads 统计信息（可选）
        preprocessing_config_file = self.main_output_dir / "01.reads_preprocessor" / "preprocessing_config.json"
        if preprocessing_config_file.exists():
            with open(preprocessing_config_file, 'r') as f:
                self.preprocessing_config = json.load(f)
        else:
            self.preprocessing_config = None
            print("  Note: Reads preprocessing config not found, some information may be unavailable")
    
    def generate_report(self):
        """
        生成可视化报告

        Returns:
            dict: 包含报告路径的字典
        """
        print(f"\n{'='*60}")
        print(f"Gap Filling Visualization")
        print(f"{'='*60}\n")

        print("[1/4] Preparing data...")
        self._prepare_data()

        print("\n[2/4] Copying sub-task reports...")
        self._copy_subtask_reports()

        print("\n[3/4] Generating HTML report...")
        html_file = self._generate_html()

        print("\n[4/4] Saving report...")
        self._save_report(html_file)

        print(f"\n{'='*60}")
        print("Visualization completed successfully!")
        print(f"{'='*60}\n")

        return {
            'output_dir': str(self.output_dir),
            'html_report': str(self.output_dir / 'Global.report.html')
        }
    
    def _prepare_data(self):
        """准备可视化数据"""
        print("  Analyzing gap filling results...")

        # 统计信息
        summary = self.gap_filling_results.get('summary', {})
        print(f"    Total gaps: {summary.get('total_gaps', 0)}")
        print(f"    Filled gaps: {summary.get('filled_gaps', 0)}")
        print(f"    Fill rate: {summary.get('fill_rate', '0%')}")

        # 按染色体组织 gaps
        gaps_by_chrom = {}
        for gap in self.gap_filling_results.get('gaps', []):
            chrom = gap['chromosome']
            if chrom not in gaps_by_chrom:
                gaps_by_chrom[chrom] = {
                    'gaps': [],
                    'total': 0,
                    'filled': 0,
                    'unfilled': 0,
                    'length': 0  # 真实染色体长度
                }
            gaps_by_chrom[chrom]['gaps'].append(gap)
            gaps_by_chrom[chrom]['total'] += 1
            if gap['filled']:
                gaps_by_chrom[chrom]['filled'] += 1
            else:
                gaps_by_chrom[chrom]['unfilled'] += 1

        # 从 integration_config 读取真实的染色体长度
        if self.integration_config:
            chrom_stats = self.integration_config.get('stats', {}).get('chromosomes', {})
            for chrom, stats in chrom_stats.items():
                if chrom in gaps_by_chrom:
                    # 使用原始长度（填充前的长度）
                    gaps_by_chrom[chrom]['length'] = stats.get('original_length', 0)

        # 如果没有 integration_config，使用 gap_end 的最大值作为估算
        for chrom, data in gaps_by_chrom.items():
            if data['length'] == 0:
                max_pos = 0
                for gap in data['gaps']:
                    gap_end = gap.get('gap_end', 0)
                    if gap_end > max_pos:
                        max_pos = gap_end
                data['length'] = max_pos

        print(f"    Chromosomes: {len(gaps_by_chrom)}")

        self.gaps_by_chrom = gaps_by_chrom

        # 计算总 BP 变化
        self.total_bp_added = 0
        if self.integration_config:
            stats = self.integration_config.get('stats', {})
            self.total_bp_added = stats.get('total_bp_added', 0)

    def _copy_subtask_reports(self):
        """复制子任务的可视化报告到总报告目录"""
        import shutil

        gap_filling_dir = self.main_output_dir / "03.gap_filling"
        copied_count = 0

        for gap in self.gap_filling_results.get('gaps', []):
            chrom = gap.get('chromosome')
            gap_id = gap.get('gap_id')
            filled = gap.get('filled', False)
            fill_direction = gap.get('fill_direction', 'none')

            # New jobs have one shared report, including unclosed results.
            if gap.get('report_direction') == 'bidirectional' or fill_direction == 'bidirectional':
                flag_to_copy = 'bidirectional'
            elif filled and fill_direction in ['left', 'right']:
                # 成功填充：复制成功的方向
                flag_to_copy = fill_direction
            else:
                # 失败：固定复制 left
                flag_to_copy = 'left'

            # 源文件路径
            source_report = gap_filling_dir / f"{chrom}.{gap_id}.{flag_to_copy}" / "visualization_output" / "gapfiller_report.html"

            if flag_to_copy == 'bidirectional':
                source_report = gap_filling_dir / f"{chrom}.{gap_id}.bidirectional" / "family-report.html"

            # 目标文件路径
            target_report = self.output_dir / f"{chrom}.{gap_id}.{flag_to_copy}.html"

            # 复制文件
            if source_report.exists():
                try:
                    if flag_to_copy == 'bidirectional':
                        # Keep source artifacts in place, but rebase report links.
                        import re
                        import os
                        from html import escape, unescape
                        from urllib.parse import quote
                        def rebase(match):
                            source = source_report.parent / unescape(match.group(1))
                            relative = os.path.relpath(source, target_report.parent)
                            return 'href="' + escape(quote(relative, safe='/.'), quote=True) + '"'
                        target_report.write_text(re.sub(r'href="([^"]+)"', rebase, source_report.read_text()))
                    else:
                        shutil.copy2(source_report, target_report)
                    copied_count += 1
                    print(f"  Copied: {chrom}.{gap_id}.{flag_to_copy}.html")
                except Exception as e:
                    print(f"  Warning: Failed to copy {source_report}: {e}")
            else:
                print(f"  Warning: Source report not found: {source_report}")

        print(f"  Total sub-task reports copied: {copied_count}")

    def _generate_html(self):
        """生成 HTML 报告"""
        print("  Creating HTML structure...")

        # 生成完整的 HTML 报告
        html_content = self._create_full_html()

        return html_content

    def _create_full_html(self):
        """创建完整的 HTML 内容"""
        summary = self.gap_filling_results.get('summary', {})

        # 准备数据
        current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        total_gaps = summary.get('total_gaps', 0)
        filled_gaps = summary.get('filled_gaps', 0)
        fill_rate = summary.get('fill_rate', '0%')
        chrom_count = len(self.gaps_by_chrom)

        bp_metric_html = ""
        if self.integration_config:
            stats = self.integration_config.get('stats', {})
            total_bp_added = getattr(self, 'total_bp_added', stats.get('total_bp_added', 0))
            bp_metric_html = f"""
            <div class="metric-card">
                <div class="metric-label">
                    <span class="lang-en">Net Length Change</span>
                    <span class="lang-zh" style="display:none;">净长度变化</span>
                </div>
                <div class="metric-value accent">{esc(self._format_bp(total_bp_added))}</div>
            </div>"""

        # 准备 JavaScript 数据
        js_summary = json.dumps(summary, indent=2)
        js_gaps_by_chrom = json.dumps(self.gaps_by_chrom, indent=2, default=str)
        js_preprocessing_config = json.dumps(self.preprocessing_config if self.preprocessing_config else {}, indent=2)
        js_integration_config = json.dumps(self.integration_config if self.integration_config else {}, indent=2)

        extra_css = """
        .table-toggle-container {
            position: sticky;
            bottom: 0;
            width: 100%;
            display: flex;
            justify-content: flex-end;
            padding: 8px 10px;
            background: rgba(255, 255, 255, 0.96);
            border-top: 1px solid var(--border);
        }

        .table-toggle-btn {
            background: var(--surface);
            color: var(--accent);
            border: 1px solid var(--accent-border);
            border-radius: var(--radius-sm);
            padding: 5px 12px;
            font-size: 0.8rem;
            font-weight: 600;
            cursor: pointer;
            transition: background-color 0.15s ease, border-color 0.15s ease;
        }

        .table-toggle-btn:hover {
            background: var(--accent-soft);
            border-color: var(--accent);
        }

        .chrom-filter-icon {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            margin-left: 6px;
            padding: 2px 5px;
            border: 1px solid transparent;
            border-radius: 4px;
            cursor: pointer;
            font-size: 0.72rem;
            color: var(--ink-muted);
            position: relative;
            user-select: none;
            vertical-align: middle;
        }

        .chrom-filter-icon:hover {
            color: var(--ink);
            background: var(--surface);
            border-color: var(--border-strong);
        }

        .chrom-filter-menu {
            display: none;
            position: absolute;
            top: calc(100% + 4px);
            left: 0;
            background: var(--surface);
            border: 1px solid var(--border-strong);
            border-radius: var(--radius-sm);
            box-shadow: 0 4px 12px rgba(15, 23, 42, 0.1);
            z-index: 1000;
            min-width: 160px;
            max-height: 280px;
            overflow-y: auto;
            text-align: left;
        }

        .chrom-filter-menu.show {
            display: block;
        }

        .chrom-filter-item {
            padding: 7px 12px;
            cursor: pointer;
            font-size: 0.84rem;
            font-weight: 400;
            color: var(--ink);
            border-bottom: 1px solid var(--border);
            white-space: nowrap;
        }

        .chrom-filter-item:last-child {
            border-bottom: none;
        }

        .chrom-filter-item:hover {
            background: var(--surface-subtle);
        }

        .chrom-filter-item.active {
            background: var(--accent-soft);
            color: var(--accent);
            font-weight: 600;
        }

        .success {
            color: var(--ok-ink);
            font-weight: 600;
        }

        .failed {
            color: var(--err-ink);
            font-weight: 600;
        }

        .chrom-selector {
            background-color: var(--surface);
            border: 1px solid var(--border-strong);
            color: var(--ink);
            padding: 6px 32px 6px 12px;
            border-radius: var(--radius-sm);
            font-size: 0.84rem;
            font-weight: 500;
            cursor: pointer;
            appearance: none;
            background-image: url('data:image/svg+xml;charset=US-ASCII,%3Csvg%20width%3D%2212%22%20height%3D%227%22%20viewBox%3D%220%200%2012%207%22%20fill%3D%22none%22%20xmlns%3D%22http%3A//www.w3.org/2000/svg%22%3E%3Cpath%20d%3D%22M1%201L6%206L11%201%22%20stroke%3D%22%23334155%22%20stroke-width%3D%221.6%22%20stroke-linecap%3D%22round%22%20stroke-linejoin%3D%22round%22/%3E%3C/svg%3E');
            background-repeat: no-repeat;
            background-position: right 10px center;
            min-width: 180px;
        }

        .chrom-selector:hover {
            border-color: var(--accent);
        }

        .track-control-bar {
            display: none;
            padding: 12px 20px;
            background: var(--surface-subtle);
            border-bottom: 1px solid var(--border);
        }

        .track-control-row {
            display: flex;
            align-items: center;
            gap: 14px;
            flex-wrap: wrap;
        }

        .track-control-label {
            font-weight: 600;
            font-size: 0.84rem;
            color: var(--ink-secondary);
            min-width: 110px;
        }

        .track-control-value {
            font-weight: 600;
            font-size: 0.85rem;
            color: var(--accent);
            font-variant-numeric: tabular-nums;
            min-width: 90px;
            text-align: right;
        }

        input[type="range"] {
            -webkit-appearance: none;
            appearance: none;
            height: 6px;
            background: var(--border-strong);
            border-radius: 3px;
            flex: 1;
            min-width: 160px;
        }

        input[type="range"]::-webkit-slider-thumb {
            -webkit-appearance: none;
            appearance: none;
            width: 16px;
            height: 16px;
            background: var(--accent);
            border: 2px solid var(--surface);
            border-radius: 50%;
            box-shadow: 0 0 0 1px var(--accent);
            cursor: pointer;
        }

        input[type="range"]::-moz-range-thumb {
            width: 16px;
            height: 16px;
            background: var(--accent);
            border: 2px solid var(--surface);
            border-radius: 50%;
            box-shadow: 0 0 0 1px var(--accent);
            cursor: pointer;
        }

        #chromVisualization {
            margin: 8px 0 12px;
            min-height: 120px;
        }

        .chrom-track {
            margin-bottom: 24px;
            padding-bottom: 18px;
            border-bottom: 1px solid var(--border);
        }

        .chrom-track:last-child {
            margin-bottom: 0;
            padding-bottom: 0;
            border-bottom: none;
        }

        .chrom-track-multirow {
            margin-bottom: 16px;
            padding-bottom: 14px;
        }

        .chrom-track-header {
            display: flex;
            justify-content: space-between;
            align-items: baseline;
            gap: 12px;
            margin-bottom: 10px;
            flex-wrap: wrap;
        }

        .chrom-track-name {
            font-family: var(--font-mono);
            font-weight: 600;
            font-size: 0.88rem;
            color: var(--ink);
            overflow-wrap: anywhere;
        }

        .chrom-track-length {
            font-size: 0.82rem;
            color: var(--ink-muted);
            font-variant-numeric: tabular-nums;
        }

        .chrom-track-bar-container {
            position: relative;
            height: 74px;
            padding: 24px 0;
        }

        .chrom-track-bar {
            position: relative;
            height: 16px;
            background: #e2e8f0;
            border: 1px solid #cbd5e1;
            border-radius: 8px;
        }

        .gap-marker {
            position: absolute;
            transform: translateX(-50%);
            cursor: pointer;
        }

        .gap-marker.top {
            top: 0;
            height: 50%;
        }

        .gap-marker.bottom {
            top: 25%;
            height: 50%;
        }

        .gap-marker-line {
            width: 2px;
            height: 100%;
            margin: 0 auto;
        }

        .gap-marker.has-width {
            transform: none;
        }

        .gap-marker.has-width .gap-marker-line {
            width: 100%;
            margin: 0;
        }

        .gap-marker.has-width .gap-marker-label {
            left: 50%;
            transform: translateX(-50%);
        }

        .gap-marker-line.filled {
            background: #15803d;
        }

        .gap-marker-line.unfilled {
            background: #b91c1c;
        }

        .gap-marker-label {
            position: absolute;
            left: 50%;
            transform: translateX(-50%);
            font-family: var(--font-mono);
            font-size: 11px;
            font-weight: 600;
            white-space: nowrap;
            cursor: pointer;
        }

        .gap-marker-label.top {
            bottom: 100%;
            margin-bottom: 3px;
        }

        .gap-marker-label.bottom {
            top: 100%;
            margin-top: 3px;
        }

        .gap-marker-label.filled {
            color: #166534;
        }

        .gap-marker-label.unfilled {
            color: #991b1b;
        }

        .gap-marker:hover .gap-marker-label {
            text-decoration: underline;
        }

        .legend {
            display: flex;
            flex-wrap: wrap;
            gap: 16px;
            margin-top: 16px;
            padding: 12px 14px;
            background: var(--surface-subtle);
            border: 1px solid var(--border);
            border-radius: var(--radius-sm);
            font-size: 0.84rem;
            color: var(--ink-secondary);
        }

        .legend-item {
            display: inline-flex;
            align-items: center;
            gap: 8px;
        }

        .legend-color {
            width: 22px;
            height: 12px;
            border-radius: 3px;
            border: 1px solid rgba(15, 23, 42, 0.15);
            flex-shrink: 0;
        }

        .tooltip {
            position: absolute;
            background: #0f172a;
            color: #f8fafc;
            border: 1px solid #334155;
            padding: 9px 12px;
            border-radius: var(--radius-sm);
            font-size: 0.78rem;
            line-height: 1.5;
            pointer-events: none;
            z-index: 1000;
            display: none;
            max-width: 280px;
            box-shadow: 0 4px 12px rgba(15, 23, 42, 0.2);
        }

        .table-link {
            font-weight: 600;
            white-space: nowrap;
        }
        """

        css = report_css(extra_css)

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AutoGapfiller - Gap Filling Report</title>
    <style>
{css}
    </style>
</head>
<body>
    <main class="container">
        <!-- Header -->
        <header class="header">
            <div class="header-top">
                <div class="header-info">
                    <div class="report-eyebrow">DEGAP · AutoGapfiller</div>
                    <h1>AutoGapfiller</h1>
                    <p class="lang-en">Gap Filling Visualization Report</p>
                    <p class="lang-zh" style="display:none;">Gap 填充可视化报告</p>
                </div>
                <div class="language-switcher" role="group" aria-label="Language Switcher">
                    <button type="button" class="language-btn active" data-lang="en" aria-pressed="true" onclick="switchLanguage('en')">EN</button>
                    <button type="button" class="language-btn" data-lang="zh" aria-pressed="false" onclick="switchLanguage('zh')">中文</button>
                </div>
            </div>
            <div class="header-meta">
                <span>
                    <span class="lang-en">Generated:</span>
                    <span class="lang-zh" style="display:none;">生成时间：</span>
                    <strong class="num">{esc(current_time)}</strong>
                </span>
            </div>
        </header>

        <!-- Summary Metrics -->
        <section class="metrics-grid" aria-label="Summary overview">
            <div class="metric-card">
                <div class="metric-label">
                    <span class="lang-en">Total Gaps</span>
                    <span class="lang-zh" style="display:none;">Gap 总数</span>
                </div>
                <div class="metric-value">{esc(total_gaps)}</div>
            </div>
            <div class="metric-card">
                <div class="metric-label">
                    <span class="lang-en">Filled Gaps</span>
                    <span class="lang-zh" style="display:none;">已填充 Gap</span>
                </div>
                <div class="metric-value ok">{esc(filled_gaps)}</div>
            </div>
            <div class="metric-card">
                <div class="metric-label">
                    <span class="lang-en">Fill Rate</span>
                    <span class="lang-zh" style="display:none;">填充率</span>
                </div>
                <div class="metric-value accent">{esc(fill_rate)}</div>
            </div>
            <div class="metric-card">
                <div class="metric-label">
                    <span class="lang-en">Chromosomes</span>
                    <span class="lang-zh" style="display:none;">涉及染色体</span>
                </div>
                <div class="metric-value">{chrom_count}</div>
            </div>{bp_metric_html}
        </section>

        <!-- Reads Information Card -->
        <section class="card">
            <div class="card-header">
                <span class="lang-en">Reads Information</span>
                <span class="lang-zh" style="display:none;">Reads 信息</span>
            </div>
            <div class="card-content">
                <div class="table-scroll">
                    <table class="kv" id="readsInfoTable">
                        <!-- Will be populated by JavaScript -->
                    </table>
                </div>
            </div>
        </section>

        <!-- Initial Genome Overview Card -->
        <section class="card">
            <div class="card-header">
                <span class="lang-en">Initial Genome Overview</span>
                <span class="lang-zh" style="display:none;">初始基因组概览</span>
            </div>
            <div class="card-content">
                <div id="initial-table-container" class="table-container-limited">
                    <table id="initialGenomeTable">
                        <thead>
                            <tr>
                                <th class="lang-en">
                                    Chromosome
                                    <span class="chrom-filter-icon" tabindex="0" role="button" aria-label="Filter chromosome" onclick="toggleFilterMenu('initial', event)" onkeydown="if(event.key==='Enter'||event.key===' '){{event.preventDefault();toggleFilterMenu('initial', event);}}">
                                        ▼
                                        <div id="initial-chrom-filter" class="chrom-filter-menu"></div>
                                    </span>
                                </th>
                                <th class="lang-zh" style="display:none;">
                                    染色体
                                    <span class="chrom-filter-icon" tabindex="0" role="button" aria-label="筛选染色体" onclick="toggleFilterMenu('initial', event)" onkeydown="if(event.key==='Enter'||event.key===' '){{event.preventDefault();toggleFilterMenu('initial', event);}}">
                                        ▼
                                        <div id="initial-chrom-filter-zh" class="chrom-filter-menu"></div>
                                    </span>
                                </th>
                                <th class="lang-en">Gap ID</th>
                                <th class="lang-zh" style="display:none;">Gap ID</th>
                                <th class="lang-en">Region</th>
                                <th class="lang-zh" style="display:none;">区域</th>
                            </tr>
                        </thead>
                        <tbody id="initialGenomeTableBody">
                            <!-- Will be populated by JavaScript -->
                        </tbody>
                    </table>
                    <div class="table-toggle-container">
                        <button type="button" id="initial-toggle-btn" class="table-toggle-btn" onclick="toggleTableExpand('initial')"><span class="lang-en">Expand</span><span class="lang-zh" style="display:none;">展开</span></button>
                    </div>
                </div>
            </div>
        </section>

        <!-- Gap Filling Tasks Card -->
        <section class="card">
            <div class="card-header">
                <span class="lang-en">Gap Filling Tasks</span>
                <span class="lang-zh" style="display:none;">Gap 填充任务</span>
            </div>
            <div class="card-content">
                <div id="tasks-table-container" class="table-container-limited">
                    <table id="gapFillingTasksTable">
                        <thead>
                            <tr>
                                <th class="lang-en">
                                    Chromosome
                                    <span class="chrom-filter-icon" tabindex="0" role="button" aria-label="Filter chromosome" onclick="toggleFilterMenu('tasks', event)" onkeydown="if(event.key==='Enter'||event.key===' '){{event.preventDefault();toggleFilterMenu('tasks', event);}}">
                                        ▼
                                        <div id="tasks-chrom-filter" class="chrom-filter-menu"></div>
                                    </span>
                                </th>
                                <th class="lang-zh" style="display:none;">
                                    染色体
                                    <span class="chrom-filter-icon" tabindex="0" role="button" aria-label="筛选染色体" onclick="toggleFilterMenu('tasks', event)" onkeydown="if(event.key==='Enter'||event.key===' '){{event.preventDefault();toggleFilterMenu('tasks', event);}}">
                                        ▼
                                        <div id="tasks-chrom-filter-zh" class="chrom-filter-menu"></div>
                                    </span>
                                </th>
                                <th class="lang-en">Gap ID</th>
                                <th class="lang-zh" style="display:none;">Gap ID</th>
                                <th class="lang-en">Filled</th>
                                <th class="lang-zh" style="display:none;">是否填补</th>
                                <th class="lang-en">Connection Type</th>
                                <th class="lang-zh" style="display:none;">连接类型</th>
                                <th class="lang-en">Direction</th>
                                <th class="lang-zh" style="display:none;">方向/模式</th>
                                <th class="lang-en">Rounds</th>
                                <th class="lang-zh" style="display:none;">轮次数量</th>
                                <th class="lang-en">Fill Length</th>
                                <th class="lang-zh" style="display:none;">填补长度</th>
                                <th class="lang-en">Details</th>
                                <th class="lang-zh" style="display:none;">详情</th>
                            </tr>
                        </thead>
                        <tbody id="gapFillingTasksTableBody">
                            <!-- Will be populated by JavaScript -->
                        </tbody>
                    </table>
                    <div class="table-toggle-container">
                        <button type="button" id="tasks-toggle-btn" class="table-toggle-btn" onclick="toggleTableExpand('tasks')"><span class="lang-en">Expand</span><span class="lang-zh" style="display:none;">展开</span></button>
                    </div>
                </div>
            </div>
        </section>

        <!-- Final Integrated Genome Card -->
        <section class="card">
            <div class="card-header">
                <span class="lang-en">Final Integrated Genome</span>
                <span class="lang-zh" style="display:none;">最终整合基因组</span>
            </div>
            <div class="card-content">
                <div id="final-table-container" class="table-container-limited">
                    <table id="finalGenomeTable">
                        <thead>
                            <tr>
                                <th class="lang-en">
                                    Chromosome
                                    <span class="chrom-filter-icon" tabindex="0" role="button" aria-label="Filter chromosome" onclick="toggleFilterMenu('final', event)" onkeydown="if(event.key==='Enter'||event.key===' '){{event.preventDefault();toggleFilterMenu('final', event);}}">
                                        ▼
                                        <div id="final-chrom-filter" class="chrom-filter-menu"></div>
                                    </span>
                                </th>
                                <th class="lang-zh" style="display:none;">
                                    染色体
                                    <span class="chrom-filter-icon" tabindex="0" role="button" aria-label="筛选染色体" onclick="toggleFilterMenu('final', event)" onkeydown="if(event.key==='Enter'||event.key===' '){{event.preventDefault();toggleFilterMenu('final', event);}}">
                                        ▼
                                        <div id="final-chrom-filter-zh" class="chrom-filter-menu"></div>
                                    </span>
                                </th>
                                <th class="lang-en">Gap ID</th>
                                <th class="lang-zh" style="display:none;">Gap ID</th>
                                <th class="lang-en">Region</th>
                                <th class="lang-zh" style="display:none;">区域</th>
                            </tr>
                        </thead>
                        <tbody id="finalGenomeTableBody">
                            <!-- Will be populated by JavaScript -->
                        </tbody>
                    </table>
                    <div class="table-toggle-container">
                        <button type="button" id="final-toggle-btn" class="table-toggle-btn" onclick="toggleTableExpand('final')"><span class="lang-en">Expand</span><span class="lang-zh" style="display:none;">展开</span></button>
                    </div>
                </div>
            </div>
        </section>

        <!-- Chromosome Visualization Card -->
        <section class="card">
            <div class="card-header">
                <span class="lang-en">Chromosome Visualization</span>
                <span class="lang-zh" style="display:none;">染色体可视化</span>
                <select id="chromSelect" class="chrom-selector" aria-label="Select chromosome" onchange="filterChromosome()">
                    <option value="all" class="lang-en">All Chromosomes</option>
                    <option value="all" class="lang-zh" style="display:none;">所有染色体</option>
                </select>
            </div>

            <!-- Track Length Control (only visible for single chromosome) -->
            <div id="trackLengthControl" class="track-control-bar">
                <div class="track-control-row">
                    <label for="trackLengthSlider" class="track-control-label">
                        <span class="lang-en">Track Length:</span>
                        <span class="lang-zh" style="display:none;">轨道长度：</span>
                    </label>
                    <input type="range" id="trackLengthSlider" min="1000000" max="300000000" step="1000000"
                           oninput="updateTrackLength()">
                    <span id="trackLengthValue" class="track-control-value">50.00 Mb</span>
                </div>
            </div>

            <div class="card-content">
                <div id="chromVisualization">
                    <!-- Will be populated by JavaScript -->
                </div>

                <div class="legend">
                    <div class="legend-item">
                        <span class="legend-color" style="background: #e2e8f0;"></span>
                        <span class="lang-en">Normal sequence</span>
                        <span class="lang-zh" style="display:none;">正常序列</span>
                    </div>
                    <div class="legend-item">
                        <span class="legend-color" style="background: #15803d;"></span>
                        <span class="lang-en">Filled gap</span>
                        <span class="lang-zh" style="display:none;">已填充 gap</span>
                    </div>
                    <div class="legend-item">
                        <span class="legend-color" style="background: #b91c1c;"></span>
                        <span class="lang-en">Unfilled gap</span>
                        <span class="lang-zh" style="display:none;">未填充 gap</span>
                    </div>
                </div>
            </div>
        </section>
    </main>

    <!-- Tooltip -->
    <div class="tooltip" id="tooltip" role="tooltip"></div>

    <script>
        // Data from Python
        const summary = {js_summary};
        const gapsByChrom = {js_gaps_by_chrom};
        const preprocessingConfig = {js_preprocessing_config};
        const integrationConfig = {js_integration_config};

        // Current language
        let currentLang = 'en';

        function escapeHtml(value) {{
            return String(value ?? '')
                .replaceAll('&', '&amp;')
                .replaceAll('<', '&lt;')
                .replaceAll('>', '&gt;')
                .replaceAll('"', '&quot;')
                .replaceAll("'", '&#039;');
        }}

        // Language switching
        function switchLanguage(lang) {{
            currentLang = lang;
            document.documentElement.setAttribute('lang', lang === 'zh' ? 'zh-CN' : 'en');

            // Update button states
            document.querySelectorAll('.language-btn').forEach(btn => {{
                const isActive = btn.getAttribute('data-lang') === lang;
                btn.classList.toggle('active', isActive);
                btn.setAttribute('aria-pressed', isActive ? 'true' : 'false');
            }});

            // Show/hide elements
            if (lang === 'en') {{
                document.querySelectorAll('.lang-en').forEach(el => el.style.display = '');
                document.querySelectorAll('.lang-zh').forEach(el => el.style.display = 'none');
            }} else {{
                document.querySelectorAll('.lang-en').forEach(el => el.style.display = 'none');
                document.querySelectorAll('.lang-zh').forEach(el => el.style.display = '');
            }}

            // Update data-i18n nodes in dynamic rows
            document.querySelectorAll('[data-i18n-en]').forEach(el => {{
                const val = el.getAttribute('data-i18n-' + lang);
                if (val !== null) {{
                    el.textContent = val;
                }}
            }});

            // Update chromosome selector default option
            const chromSelect = document.getElementById('chromSelect');
            if (chromSelect) {{
                const allOpt = chromSelect.querySelector('option[value="all"]');
                if (allOpt) {{
                    allOpt.textContent = lang === 'zh' ? '所有染色体' : 'All Chromosomes';
                }}
            }}
        }}

        // Expand/Collapse table
        function toggleTableExpand(which) {{
            const containerMap = {{
                'initial': 'initial-table-container',
                'tasks': 'tasks-table-container',
                'final': 'final-table-container'
            }};

            const container = document.getElementById(containerMap[which]);
            const btn = document.getElementById(which + '-toggle-btn');
            const expanded = container.classList.toggle('expanded');

            // Update button text based on current language
            if (currentLang === 'en') {{
                btn.innerHTML = expanded ? '<span class="lang-en">Collapse</span><span class="lang-zh" style="display:none;">收起</span>' : '<span class="lang-en">Expand</span><span class="lang-zh" style="display:none;">展开</span>';
            }} else {{
                btn.innerHTML = expanded ? '<span class="lang-en" style="display:none;">Collapse</span><span class="lang-zh">收起</span>' : '<span class="lang-en" style="display:none;">Expand</span><span class="lang-zh">展开</span>';
            }}
        }}

        // Format number with commas
        function formatNumber(num) {{
            return String(num).replace(/\\B(?=(\\d{{3}})+(?!\\d))/g, ",");
        }}

        // Format BP
        function formatBP(bp) {{
            if (Math.abs(bp) >= 1000000) {{
                return (bp / 1000000).toFixed(2) + ' Mb';
            }} else if (Math.abs(bp) >= 1000) {{
                return (bp / 1000).toFixed(2) + ' Kb';
            }} else {{
                return bp + ' bp';
            }}
        }}

        // Populate Reads Information
        function populateReadsInfo() {{
            const table = document.getElementById('readsInfoTable');

            if (!preprocessingConfig || Object.keys(preprocessingConfig).length === 0) {{
                table.innerHTML = '<tr><td colspan="2" class="empty-state"><span class="lang-en">No reads information available</span><span class="lang-zh" style="display:none;">暂无 Reads 信息</span></td></tr>';
                return;
            }}

            const dataType = String(preprocessingConfig.data_type || 'unknown');
            const hifiStats = preprocessingConfig.hifi_stats || {{}};
            const ontStats = preprocessingConfig.ont_stats || {{}};

            let html = '';

            // Data Type
            const displayType = dataType.charAt(0).toUpperCase() + dataType.slice(1);
            html += `
                <tr>
                    <td><span class="lang-en">Data Type</span><span class="lang-zh" style="display:none;">数据类型</span></td>
                    <td>${{escapeHtml(displayType)}}</td>
                </tr>
            `;

            // HiFi Stats
            if (hifiStats && Object.keys(hifiStats).length > 0) {{
                html += `
                    <tr>
                        <td>HiFi Reads</td>
                        <td class="num">${{escapeHtml(formatNumber(hifiStats.Number || hifiStats.NumSeqs || 0))}} reads</td>
                    </tr>
                    <tr>
                        <td>HiFi Max Length</td>
                        <td class="num">${{escapeHtml(formatNumber(hifiStats.MaxLength || 0))}} bp</td>
                    </tr>
                `;
            }}

            // ONT Stats
            if (ontStats && Object.keys(ontStats).length > 0) {{
                html += `
                    <tr>
                        <td>ONT Reads</td>
                        <td class="num">${{escapeHtml(formatNumber(ontStats.Number || ontStats.NumSeqs || 0))}} reads</td>
                    </tr>
                    <tr>
                        <td>ONT Max Length</td>
                        <td class="num">${{escapeHtml(formatNumber(ontStats.MaxLength || 0))}} bp</td>
                    </tr>
                `;
            }}

            // Seed Length - use unified seed_length (max of HiFi and ONT)
            const seedLength = preprocessingConfig.seed_length || 0;
            if (seedLength > 0) {{
                html += `
                    <tr>
                        <td><span class="lang-en">Seed Length</span><span class="lang-zh" style="display:none;">Seed 长度</span></td>
                        <td class="num">${{escapeHtml(formatNumber(seedLength))}} bp</td>
                    </tr>
                `;
            }}

            table.innerHTML = html;
        }}

        // Populate Initial Genome Table
        function populateInitialGenomeTable() {{
            const tbody = document.getElementById('initialGenomeTableBody');

            let html = '';
            let rowCount = 0;

            // Sort chromosomes
            const chromNames = Object.keys(gapsByChrom).sort();

            chromNames.forEach(chrom => {{
                const chromData = gapsByChrom[chrom];

                // Sort gaps by gap_id
                const sortedGaps = chromData.gaps.sort((a, b) => a.gap_id - b.gap_id);

                sortedGaps.forEach(gap => {{
                    const region = `${{formatNumber(gap.gap_start)}} - ${{formatNumber(gap.gap_end)}}`;

                    html += `
                        <tr>
                            <td class="mono">${{escapeHtml(chrom)}}</td>
                            <td class="num">${{escapeHtml(gap.gap_id)}}</td>
                            <td class="num">${{escapeHtml(region)}}</td>
                        </tr>
                    `;
                    rowCount++;
                }});
            }});

            if (rowCount === 0) {{
                html = '<tr><td colspan="3" class="empty-state"><span class="lang-en">No initial gaps recorded</span><span class="lang-zh" style="display:none;">暂无初始 Gap 记录</span></td></tr>';
            }}

            tbody.innerHTML = html;

            // Show/hide toggle button based on row count (show if >= 6 rows)
            const toggleContainer = document.querySelector('#initial-table-container .table-toggle-container');
            if (rowCount >= 6) {{
                toggleContainer.style.display = 'flex';
            }} else {{
                toggleContainer.style.display = 'none';
            }}
        }}

        // Populate Gap Filling Tasks Table
        function populateGapFillingTasksTable() {{
            const tbody = document.getElementById('gapFillingTasksTableBody');

            let html = '';
            let rowCount = 0;

            // Sort chromosomes
            const chromNames = Object.keys(gapsByChrom).sort();

            chromNames.forEach(chrom => {{
                const chromData = gapsByChrom[chrom];

                // Sort gaps by gap_id
                const sortedGaps = chromData.gaps.sort((a, b) => a.gap_id - b.gap_id);

                sortedGaps.forEach(gap => {{
                    // Extended - use badge style
                    const extended = gap.filled ?
                        '<span class="badge badge-yes" data-i18n-en="Yes" data-i18n-zh="是">Yes</span>' :
                        '<span class="badge badge-no" data-i18n-en="No" data-i18n-zh="否">No</span>';

                    // Connection Type - use badge style
                    let connectionType = '<span class="badge badge-na">-</span>';
                    if (gap.filled) {{
                        if (gap.fill_type === 'extension') {{
                            connectionType = '<span class="badge badge-extension" data-i18n-en="Extension" data-i18n-zh="延伸">Extension</span>';
                        }} else if (gap.fill_type === 'direct_connection') {{
                            connectionType = '<span class="badge badge-direct" data-i18n-en="Direct" data-i18n-zh="直连">Direct</span>';
                        }}
                    }}

                    // Flag column - show fill_direction (left/right/none)
                    const displayFlag = gap.fill_direction || 'none';

                    const rounds = gap.filled ? gap.rounds : '-';

                    // Fill length - show formatted BP if filled
                    const hasFillLen = Boolean(gap.filled && gap.fill_length);
                    const fillLength = hasFillLen ? formatBP(gap.fill_length) : '-';
                    const fillLenClass = hasFillLen ? 'num success' : 'num muted';

                    // Details - View link to sub-task report
                    // For filled gaps: use fill_direction, for unfilled: use 'left'
                    const linkFlag = gap.report_direction || (gap.filled && gap.fill_direction ? gap.fill_direction : 'left');
                    const reportFile = `${{chrom}}.${{gap.gap_id}}.${{linkFlag}}.html`;
                    const details = `<a class="table-link" href="${{escapeHtml(reportFile)}}" target="_blank" rel="noopener" title="View details for ${{escapeHtml(chrom)}} Gap ${{escapeHtml(gap.gap_id)}}" data-i18n-en="View" data-i18n-zh="查看">View</a>`;

                    html += `
                        <tr>
                            <td class="mono">${{escapeHtml(chrom)}}</td>
                            <td class="num">${{escapeHtml(gap.gap_id)}}</td>
                            <td>${{extended}}</td>
                            <td>${{connectionType}}</td>
                            <td class="mono">${{escapeHtml(displayFlag)}}</td>
                            <td class="num">${{escapeHtml(rounds)}}</td>
                            <td class="${{fillLenClass}}">${{escapeHtml(fillLength)}}</td>
                            <td>${{details}}</td>
                        </tr>
                    `;
                    rowCount++;
                }});
            }});

            if (rowCount === 0) {{
                html = '<tr><td colspan="8" class="empty-state"><span class="lang-en">No gap filling tasks recorded</span><span class="lang-zh" style="display:none;">暂无 Gap 填充任务记录</span></td></tr>';
            }}

            tbody.innerHTML = html;

            // Show/hide toggle button based on row count (show if >= 6 rows)
            const toggleContainer = document.querySelector('#tasks-table-container .table-toggle-container');
            if (rowCount >= 6) {{
                toggleContainer.style.display = 'flex';
            }} else {{
                toggleContainer.style.display = 'none';
            }}
        }}

        // Populate Final Integrated Genome Table
        function populateFinalGenomeTable() {{
            const tbody = document.getElementById('finalGenomeTableBody');

            let html = '';
            let rowCount = 0;

            // Sort chromosomes
            const chromNames = Object.keys(gapsByChrom).sort();

            chromNames.forEach(chrom => {{
                const chromData = gapsByChrom[chrom];

                // Sort all gaps by gap_id
                const allGaps = chromData.gaps.sort((a, b) => a.gap_id - b.gap_id);

                // Calculate position offset based on filled gaps
                let positionOffset = 0;

                allGaps.forEach(gap => {{
                    if (!gap.filled) {{
                        // This gap still exists in the final genome
                        // Calculate new position by adding offset from all previous filled gaps
                        const newStart = gap.gap_start + positionOffset;
                        const newEnd = gap.gap_end + positionOffset;
                        const region = `${{formatNumber(newStart)}} - ${{formatNumber(newEnd)}}`;

                        html += `
                            <tr>
                                <td class="mono">${{escapeHtml(chrom)}}</td>
                                <td class="num">${{escapeHtml(gap.gap_id)}}</td>
                                <td class="num">${{escapeHtml(region)}}</td>
                            </tr>
                        `;
                        rowCount++;
                    }} else {{
                        // This gap was filled, update offset for subsequent gaps
                        // Offset = fill_length - original_gap_length
                        const originalGapLen = gap.gap_end - gap.gap_start;
                        const fillLength = gap.fill_length || 0;
                        positionOffset += (fillLength - originalGapLen);
                    }}
                }});
            }});

            // If no unfilled gaps, show a message
            if (html === '') {{
                html = '<tr><td colspan="3" class="empty-state success"><span data-i18n-en="All gaps have been filled." data-i18n-zh="所有 gaps 已填充。">All gaps have been filled.</span></td></tr>';
            }}

            tbody.innerHTML = html;

            // Show/hide toggle button based on row count (show if >= 6 rows)
            const toggleContainer = document.querySelector('#final-table-container .table-toggle-container');
            if (rowCount >= 6) {{
                toggleContainer.style.display = 'flex';
            }} else {{
                toggleContainer.style.display = 'none';
            }}
        }}

        // Populate Chromosome Selector
        function populateChromSelector() {{
            const select = document.getElementById('chromSelect');
            let selectHtml = `<option value="all">${{currentLang === 'zh' ? '所有染色体' : 'All Chromosomes'}}</option>`;

            // Sort chromosomes
            const chromNames = Object.keys(gapsByChrom).sort();

            chromNames.forEach(chrom => {{
                selectHtml += `<option value="${{escapeHtml(chrom)}}">${{escapeHtml(chrom)}}</option>`;
            }});

            select.innerHTML = selectHtml;
        }}

        // Current track length for multi-row view
        let currentTrackLength = 0;

        // Draw Chromosome Visualization
        function drawChromVisualization(selectedChrom = 'all') {{
            const container = document.getElementById('chromVisualization');
            container.innerHTML = '';

            const chromNames = selectedChrom === 'all'
                ? Object.keys(gapsByChrom).sort()
                : [selectedChrom];

            if (chromNames.length === 0) {{
                container.innerHTML = '<div class="empty-state"><span class="lang-en">No chromosomes to display</span><span class="lang-zh" style="display:none;">暂无可显示的染色体</span></div>';
                return;
            }}

            // Show/hide track length control
            const trackLengthControl = document.getElementById('trackLengthControl');
            if (selectedChrom !== 'all' && gapsByChrom[selectedChrom]) {{
                trackLengthControl.style.display = 'block';
                const chromData = gapsByChrom[selectedChrom];
                const chromLength = chromData.length;

                // Set slider range and default value
                const slider = document.getElementById('trackLengthSlider');
                slider.max = chromLength;
                slider.min = Math.min(1000000, chromLength); // 1 Mb minimum

                // Default: chromosome length / 5
                if (currentTrackLength === 0 || currentTrackLength > chromLength) {{
                    currentTrackLength = Math.floor(chromLength / 5);
                }}
                slider.value = currentTrackLength;
                updateTrackLengthDisplay(currentTrackLength);
            }} else {{
                trackLengthControl.style.display = 'none';
                currentTrackLength = 0;
            }}

            // Find max chromosome length for scaling
            let maxChromLength = 0;
            chromNames.forEach(chrom => {{
                const chromData = gapsByChrom[chrom];
                if (!chromData) return;
                const chromLength = chromData.length;
                if (chromLength > maxChromLength) {{
                    maxChromLength = chromLength;
                }}
            }});

            chromNames.forEach(chrom => {{
                const chromData = gapsByChrom[chrom];
                if (!chromData) return;
                const chromLength = chromData.length;

                // Multi-row view for single chromosome
                if (selectedChrom !== 'all' && currentTrackLength > 0) {{
                    drawMultiRowChromosome(chrom, chromData, chromLength, currentTrackLength, container);
                }} else {{
                    // Single row view for all chromosomes
                    drawSingleRowChromosome(chrom, chromData, chromLength, maxChromLength || 1, container);
                }}
            }});
        }}

        // Draw single row chromosome (for "All Chromosomes" view)
        function drawSingleRowChromosome(chrom, chromData, chromLength, maxChromLength, container) {{
            // Create chromosome track
            const trackDiv = document.createElement('div');
            trackDiv.className = 'chrom-track';

            // Header
            const headerDiv = document.createElement('div');
            headerDiv.className = 'chrom-track-header';
            headerDiv.innerHTML = `
                <span class="chrom-track-name">${{escapeHtml(chrom)}}</span>
                <span class="chrom-track-length">${{escapeHtml(formatBP(chromLength))}}</span>
            `;
            trackDiv.appendChild(headerDiv);

            // Bar container
            const barContainerDiv = document.createElement('div');
            barContainerDiv.className = 'chrom-track-bar-container';

            // Bar (scaled relative to max chromosome)
            const barDiv = document.createElement('div');
            barDiv.className = 'chrom-track-bar';
            const barWidth = (chromLength / maxChromLength * 100).toFixed(2);
            barDiv.style.width = barWidth + '%';

            barContainerDiv.appendChild(barDiv);

            // Calculate bar width percentage for gap positioning
            const barWidthPercent = (chromLength / maxChromLength * 100);

            // Add gap markers
            chromData.gaps.forEach((gap, index) => {{
                addGapMarker(barContainerDiv, gap, index, chromLength || 1, chrom, null, false, barWidthPercent);
            }});

            trackDiv.appendChild(barContainerDiv);
            container.appendChild(trackDiv);
        }}

        // Draw multi-row chromosome (for single chromosome view)
        function drawMultiRowChromosome(chrom, chromData, chromLength, trackLength, container) {{
            const numRows = Math.ceil(chromLength / trackLength);

            for (let row = 0; row < numRows; row++) {{
                const rowStart = row * trackLength;
                const rowEnd = Math.min((row + 1) * trackLength, chromLength);
                const rowLength = rowEnd - rowStart;

                // Create chromosome track
                const trackDiv = document.createElement('div');
                trackDiv.className = 'chrom-track chrom-track-multirow';

                // Header
                const headerDiv = document.createElement('div');
                headerDiv.className = 'chrom-track-header';
                headerDiv.innerHTML = `
                    <span class="chrom-track-name">${{escapeHtml(chrom)}} (Row ${{row + 1}}/${{numRows}})</span>
                    <span class="chrom-track-length">${{escapeHtml(formatBP(rowStart))}} - ${{escapeHtml(formatBP(rowEnd))}}</span>
                `;
                trackDiv.appendChild(headerDiv);

                // Bar container
                const barContainerDiv = document.createElement('div');
                barContainerDiv.className = 'chrom-track-bar-container';

                // Bar (width based on actual length, not full track length for last row)
                const barDiv = document.createElement('div');
                barDiv.className = 'chrom-track-bar';
                const barWidth = (rowLength / trackLength * 100).toFixed(2);
                barDiv.style.width = barWidth + '%';

                barContainerDiv.appendChild(barDiv);

                // Add gap markers for gaps in this row
                chromData.gaps.forEach((gap, index) => {{
                    if (gap.gap_start >= rowStart && gap.gap_start < rowEnd) {{
                        // Calculate position relative to this row
                        const relativePosition = gap.gap_start - rowStart;
                        addGapMarker(barContainerDiv, gap, index, trackLength, chrom, relativePosition, true);
                    }}
                }});

                trackDiv.appendChild(barContainerDiv);
                container.appendChild(trackDiv);
            }}
        }}

        // Add gap marker to a bar container
        function addGapMarker(barContainerDiv, gap, index, chromLength, chrom, customPosition = null, useVariableWidth = false, barWidthPercent = 100) {{
            const isTop = index % 2 === 0;

            const markerDiv = document.createElement('div');
            markerDiv.className = `gap-marker ${{isTop ? 'top' : 'bottom'}}`;
            markerDiv.setAttribute('title', `${{chrom}} Gap ${{gap.gap_id}}: ${{gap.filled ? 'Filled' : 'Unfilled'}}`);

            // Position based on gap_start (or custom position for multi-row)
            const gapPosition = customPosition !== null ? customPosition : gap.gap_start;

            if (useVariableWidth) {{
                // Multi-row view: draw gap as a region with width
                // Use fill_length for filled gaps, gap_len for unfilled gaps
                let gapDisplayLength;
                if (gap.filled && gap.fill_length > 0) {{
                    // Filled gap: use fill_length (extended length)
                    gapDisplayLength = gap.fill_length;
                }} else {{
                    // Unfilled gap or direct connection (negative): use original gap_len
                    gapDisplayLength = gap.gap_len;
                }}

                const gapWidthPercent = (gapDisplayLength / chromLength * 100);
                const minWidthPercent = 0.1; // Minimum 0.1% width to ensure visibility
                const finalWidthPercent = Math.max(gapWidthPercent, minWidthPercent);

                markerDiv.classList.add('has-width');
                markerDiv.style.left = (gapPosition / chromLength * 100).toFixed(2) + '%';
                markerDiv.style.width = finalWidthPercent.toFixed(2) + '%';
            }} else {{
                // Single-row view: draw gap as a line
                // Position gap relative to the bar's actual width, not the full container
                const gapPositionInBar = (gapPosition / chromLength * 100);
                const finalPosition = gapPositionInBar * barWidthPercent / 100;
                markerDiv.style.left = finalPosition.toFixed(2) + '%';
            }}

            // Line (extends from track to label direction)
            const lineDiv = document.createElement('div');
            lineDiv.className = `gap-marker-line ${{gap.filled ? 'filled' : 'unfilled'}}`;
            markerDiv.appendChild(lineDiv);

            // Label (alternate top/bottom, outside the track)
            const labelDiv = document.createElement('div');
            labelDiv.className = `gap-marker-label ${{isTop ? 'top' : 'bottom'}} ${{gap.filled ? 'filled' : 'unfilled'}}`;
            labelDiv.textContent = `${{gap.gap_id}}`;
            markerDiv.appendChild(labelDiv);

            // Tooltip on hover
            markerDiv.addEventListener('mouseenter', (e) => {{
                showTooltip(e, gap, chrom);
            }});
            markerDiv.addEventListener('mouseleave', hideTooltip);

            barContainerDiv.appendChild(markerDiv);
        }}

        // Update track length display
        function updateTrackLengthDisplay(length) {{
            const valueSpan = document.getElementById('trackLengthValue');
            valueSpan.textContent = formatBP(length);
        }}

        // Update track length from slider
        function updateTrackLength() {{
            const slider = document.getElementById('trackLengthSlider');
            currentTrackLength = parseInt(slider.value);
            updateTrackLengthDisplay(currentTrackLength);

            // Redraw visualization
            const select = document.getElementById('chromSelect');
            const selectedChrom = select.value;
            drawChromVisualization(selectedChrom);
        }}

        // Show tooltip
        function showTooltip(event, gap, chrom) {{
            const tooltip = document.getElementById('tooltip');
            const status = gap.filled ? 'Filled' : 'Unfilled';
            const statusZh = gap.filled ? '已填充' : '未填充';

            let content = `
                <strong>${{escapeHtml(chrom)}}.gap${{escapeHtml(gap.gap_id)}}</strong><br>
                <span class="lang-en">Status: ${{status}}</span>
                <span class="lang-zh" style="display:none;">状态：${{statusZh}}</span><br>
                <span class="lang-en">Position: ${{escapeHtml(formatNumber(gap.gap_start))}} - ${{escapeHtml(formatNumber(gap.gap_end))}}</span>
                <span class="lang-zh" style="display:none;">位置：${{escapeHtml(formatNumber(gap.gap_start))}} - ${{escapeHtml(formatNumber(gap.gap_end))}}</span><br>
                <span class="lang-en">Initial Length: ${{escapeHtml(formatBP(gap.gap_len))}}</span>
                <span class="lang-zh" style="display:none;">初始长度：${{escapeHtml(formatBP(gap.gap_len))}}</span><br>
            `;

            // Extended length (only for filled gaps)
            if (gap.filled && gap.fill_length > 0) {{
                content += `
                    <span class="lang-en">Extended Length: ${{escapeHtml(formatBP(gap.fill_length))}}</span>
                    <span class="lang-zh" style="display:none;">延伸后长度：${{escapeHtml(formatBP(gap.fill_length))}}</span><br>
                `;
            }} else {{
                content += `
                    <span class="lang-en">Extended Length: -</span>
                    <span class="lang-zh" style="display:none;">延伸后长度：-</span><br>
                `;
            }}

            if (gap.filled) {{
                content += `<span class="lang-en">Type: ${{escapeHtml(gap.fill_type)}}</span>`;
                content += `<span class="lang-zh" style="display:none;">类型：${{escapeHtml(gap.fill_type)}}</span>`;
            }}

            tooltip.innerHTML = content;
            tooltip.style.display = 'block';
            tooltip.style.left = event.pageX + 10 + 'px';
            tooltip.style.top = event.pageY + 10 + 'px';

            // Update language display
            if (currentLang === 'zh') {{
                tooltip.querySelectorAll('.lang-en').forEach(el => el.style.display = 'none');
                tooltip.querySelectorAll('.lang-zh').forEach(el => el.style.display = '');
            }}
        }}

        // Hide tooltip
        function hideTooltip() {{
            document.getElementById('tooltip').style.display = 'none';
        }}

        // Filter chromosome
        function filterChromosome() {{
            const select = document.getElementById('chromSelect');
            const selectedChrom = select.value;
            drawChromVisualization(selectedChrom);
        }}

        // Chromosome filter functions
        let currentFilterTable = null;
        const selectedChroms = {{
            initial: 'all',
            tasks: 'all',
            final: 'all'
        }};

        function buildFilterMenu(container, tableType, allLabel, chromNames) {{
            container.innerHTML = '';
            const allItem = document.createElement('div');
            allItem.className = 'chrom-filter-item active';
            allItem.textContent = allLabel;
            allItem.addEventListener('click', (e) => {{
                e.stopPropagation();
                selectChrom(tableType, 'all');
            }});
            container.appendChild(allItem);

            chromNames.forEach(chrom => {{
                const item = document.createElement('div');
                item.className = 'chrom-filter-item';
                item.textContent = chrom;
                item.addEventListener('click', (e) => {{
                    e.stopPropagation();
                    selectChrom(tableType, chrom);
                }});
                container.appendChild(item);
            }});
        }}

        function initChromFilters() {{
            const chromNames = Object.keys(gapsByChrom).sort();

            // Initialize filters for all three tables
            ['initial', 'tasks', 'final'].forEach(tableType => {{
                const containerEn = document.getElementById(`${{tableType}}-chrom-filter`);
                const containerZh = document.getElementById(`${{tableType}}-chrom-filter-zh`);
                buildFilterMenu(containerEn, tableType, 'All', chromNames);
                buildFilterMenu(containerZh, tableType, '全部', chromNames);
            }});
        }}

        function toggleFilterMenu(tableType, event) {{
            event.stopPropagation();

            // Close all other menus
            document.querySelectorAll('.chrom-filter-menu').forEach(menu => {{
                if (!menu.id.includes(tableType)) {{
                    menu.classList.remove('show');
                }}
            }});

            // Toggle current menu
            const menuEn = document.getElementById(`${{tableType}}-chrom-filter`);
            const menuZh = document.getElementById(`${{tableType}}-chrom-filter-zh`);

            menuEn.classList.toggle('show');
            menuZh.classList.toggle('show');

            currentFilterTable = menuEn.classList.contains('show') ? tableType : null;
        }}

        function selectChrom(tableType, selectedChrom) {{
            selectedChroms[tableType] = selectedChrom;

            // Update active state
            const menuEn = document.getElementById(`${{tableType}}-chrom-filter`);
            const menuZh = document.getElementById(`${{tableType}}-chrom-filter-zh`);

            [menuEn, menuZh].forEach(menu => {{
                menu.querySelectorAll('.chrom-filter-item').forEach(item => {{
                    item.classList.remove('active');
                }});
            }});

            // Set active on selected items
            const items = document.querySelectorAll(`#${{tableType}}-chrom-filter .chrom-filter-item, #${{tableType}}-chrom-filter-zh .chrom-filter-item`);
            items.forEach(item => {{
                const text = item.textContent.trim();
                if ((selectedChrom === 'all' && (text === 'All' || text === '全部')) || text === selectedChrom) {{
                    item.classList.add('active');
                }}
            }});

            // Filter table
            filterTableByChrom(tableType, selectedChrom);

            // Close menu
            menuEn.classList.remove('show');
            menuZh.classList.remove('show');
        }}

        function filterTableByChrom(tableType, selectedChrom) {{
            const tableId = tableType === 'initial' ? 'initialGenomeTable' :
                           tableType === 'tasks' ? 'gapFillingTasksTable' :
                           'finalGenomeTable';

            const tbody = document.getElementById(tableId + 'Body');
            const rows = tbody.querySelectorAll('tr');

            rows.forEach(row => {{
                const chromCell = row.cells[0];
                if (chromCell) {{
                    const chromName = chromCell.textContent.trim();
                    if (selectedChrom === 'all' || chromName === selectedChrom) {{
                        row.style.display = '';
                    }} else {{
                        row.style.display = 'none';
                    }}
                }}
            }});
        }}

        // Close filter menu when clicking outside
        document.addEventListener('click', (event) => {{
            if (!event.target.closest('.chrom-filter-icon')) {{
                document.querySelectorAll('.chrom-filter-menu').forEach(menu => {{
                    menu.classList.remove('show');
                }});
                currentFilterTable = null;
            }}
        }});

        // Initialize
        document.addEventListener('DOMContentLoaded', () => {{
            populateReadsInfo();
            populateInitialGenomeTable();
            populateGapFillingTasksTable();
            populateFinalGenomeTable();
            populateChromSelector();
            drawChromVisualization();
            initChromFilters();
        }});
    </script>
</body>
</html>
"""

        return html

    def _format_bp(self, bp):
        """格式化碱基对数量"""
        if abs(bp) >= 1000000:
            return f"{bp/1000000:+.2f} Mb"
        elif abs(bp) >= 1000:
            return f"{bp/1000:+.2f} Kb"
        else:
            return f"{bp:+,} bp"

    def _save_report(self, html_content):
        """保存 HTML 报告"""
        html_file = self.output_dir / 'Global.report.html'

        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        print(f"  Report saved: {html_file}")


def main():
    """命令行接口"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate gap filling visualization report',
        epilog='''
Example:
  %(prog)s -o output

Note: This script requires gap filling to be completed first.
      It will read results from:
      - output/03.gap_filling/gap_filling_results.json (required)
      - output/04.genome_integration/integration_config.json (optional)
      - output/01.reads_preprocessor/preprocessing_config.json (optional)

Features:
  - Interactive HTML report with chromosome visualization
  - Bilingual support (English/Chinese)
  - Gap filling statistics and distribution
  - Chromosome band visualization with gap positions
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--output', '-o', required=True, help='Output directory')

    args = parser.parse_args()

    try:
        visualizer = GapFillingVisualizer(args.output)
        result = visualizer.generate_report()

        print("\nVisualization result:")
        print(json.dumps(result, indent=2, default=str))
        print(f"\n✅ Open the HTML file in a web browser to view the report:")
        print(f"   {result['html_report']}")
        print(f"\n💡 The report supports:")
        print(f"   - Language switching (EN/中文)")
        print(f"   - Interactive chromosome visualization")
        print(f"   - Chromosome filtering (dropdown menu)")
        print(f"   - Gap details on hover")

    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()

