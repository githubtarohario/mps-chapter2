# -*- coding: utf-8 -*-
"""Markdown -> 印刷用 HTML (日本語・表・罫線アート対応)"""
import sys, io, re, markdown


def sub_sup(s):
    """数式テキストを HTML の上付き・下付きに変換する。

    markdown が出力した時点で HTML エスケープ済みなので、
    そのまま正規表現をかけてよい。
      _{xxx} -> <sub>xxx</sub>      ^{xxx} -> <sup>xxx</sup>
      _x     -> <sub>x</sub>        ^x     -> <sup>x</sup>
      ⁰ ¹ ² ³ も <sup> に統一して字体を揃える
    """
    s = re.sub(r"_\{([^}]*)\}", r"<sub>\1</sub>", s)
    s = re.sub(r"\^\{([^}]*)\}", r"<sup>\1</sup>", s)
    s = re.sub(r"_([A-Za-z0-9])", r"<sub>\1</sub>", s)
    s = re.sub(r"\^([A-Za-z0-9*])", r"<sup>\1</sup>", s)
    for ch, n in (("⁰", "0"), ("¹", "1"), ("²", "2"), ("³", "3")):
        s = s.replace(ch, "<sup>%s</sup>" % n)
    return s


def render_formula(m):
    """```formula ブロック -> 数式用の div

    行末に「…注釈」があれば、注釈だけゴシック体の小さい字で組む。
    """
    rows = []
    for line in m.group(1).split("\n"):
        if not line.strip():
            continue
        if "…" in line:                       # …
            fx, note = line.split("…", 1)
            rows.append('<div class="fline"><span class="fx">%s</span>'
                        '<span class="fnote">%s</span></div>'
                        % (sub_sup(fx.strip()), note.strip()))
        else:
            rows.append('<div class="fline"><span class="fx">%s</span></div>'
                        % sub_sup(line.strip()))
    return '<div class="formula">%s</div>' % "".join(rows)


def render_cases(m):
    """```cases ブロック -> 場合分け（大きな中括弧つき）

    1 行目が左辺、2 行目以降が「値 @条件」。
    """
    lines = [l for l in m.group(1).split("\n") if l.strip()]
    lhs = sub_sup(lines[0].strip())
    rows = []
    for line in lines[1:]:
        if "@" in line:
            val, cond = line.split("@", 1)
            rows.append('<div><span class="cv">%s</span><span class="cc">%s</span></div>'
                        % (sub_sup(val.strip()), sub_sup(cond.strip())))
        else:
            rows.append('<div><span class="cv">%s</span></div>' % sub_sup(line.strip()))
    return ('<div class="formula cases"><span class="fx">%s</span>'
            '<span class="brace">{</span><span class="rows">%s</span></div>'
            % (lhs, "".join(rows)))

src, dst, title = sys.argv[1], sys.argv[2], sys.argv[3]

with io.open(src, encoding="utf-8-sig") as f:
    text = f.read()

html_body = markdown.markdown(
    text,
    # nl2br は使わない。原文は読みやすさのために手で折り返してあるが、
    # PDF では紙幅に合わせて再流動させたほうがきれいに組める。
    extensions=["tables", "fenced_code", "sane_lists", "attr_list"],
    output_format="html5",
)

# ```formula / ```cases を数式用の組版に差し替える
html_body = re.sub(r'<pre><code class="language-formula">(.*?)</code></pre>',
                   render_formula, html_body, flags=re.S)
html_body = re.sub(r'<pre><code class="language-cases">(.*?)</code></pre>',
                   render_cases, html_body, flags=re.S)

# 表の中などに書いた <span class="mi">…</span> は行中の数式として組む
html_body = re.sub(r'<span class="mi">(.*?)</span>',
                   lambda m: '<span class="mi">%s</span>' % sub_sup(m.group(1)),
                   html_body, flags=re.S)

CSS = """
@page { size: A4; margin: 15mm 14mm 18mm 14mm; }

html { -webkit-print-color-adjust: exact; print-color-adjust: exact; }

body {
  font-family: "Yu Gothic","Meiryo","MS PGothic","Segoe UI Symbol",sans-serif;
  font-size: 10pt;
  line-height: 1.75;
  color: #1a1a1a;
  margin: 0;
}

/* ---- 見出し ---- */
h1 {
  font-size: 20pt; line-height: 1.4; margin: 0 0 6px;
  padding-bottom: 8px; border-bottom: 3px solid #2f5d8a; color: #17334f;
}
h2 {
  font-size: 15pt; margin: 0 0 12px; padding: 6px 0 6px 10px;
  border-left: 6px solid #2f5d8a; background: #eef3f8; color: #17334f;
  break-before: page; break-after: avoid;
}
h2:first-of-type { break-before: avoid; }
h3 {
  font-size: 12.5pt; margin: 18px 0 8px; padding-bottom: 3px;
  border-bottom: 1px solid #c3cedb; color: #1d4266; break-after: avoid;
}
h4 {
  font-size: 11pt; margin: 14px 0 6px; color: #2a4a6a; break-after: avoid;
}

p { margin: 7px 0; }
strong { color: #0f2b45; }

/* ---- コードブロック（罫線アートの桁合わせが要なので等幅CJKフォント） ---- */
pre {
  font-family: "MS Gothic","Osaka-Mono",monospace;
  font-size: 8.5pt; line-height: 1.35;
  background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 4px;
  padding: 8px 10px; margin: 9px 0;
  white-space: pre; overflow: hidden;
  /* 途中でページが変わらないようにする。1 ページに収まらない長い
     ブロックはブラウザが自動的に分割するので、図が壊れる心配はない */
  break-inside: avoid;
}
pre code { font-family: inherit; font-size: inherit; background: none; padding: 0; border: none; }
code {
  font-family: "MS Gothic",monospace; font-size: 9pt;
  background: #eef1f4; border: 1px solid #dde2e7; border-radius: 3px;
  padding: 0 3px;
}

/* ---- 数式（Cambria Math で組む。日本語は游ゴシックにフォールバック） ---- */
.formula {
  font-family: "Cambria Math","Cambria","Times New Roman","Yu Gothic",serif;
  font-size: 12pt; line-height: 1.95;
  margin: 11px 0 11px 6px; padding: 5px 0 5px 15px;
  border-left: 3px solid #b9c7d6;
  break-inside: avoid;
}
.formula .fline { white-space: normal; }
.formula .fx { display: inline-block; min-width: 56%; }
.formula .fnote {
  font-family: "Yu Gothic","Meiryo",sans-serif;
  font-size: 9pt; color: #55636f;
}
/* 上付き・下付きが行間を広げないようにする */
.formula sub, .formula sup { font-size: 66%; line-height: 0; position: relative; }
.formula sub { bottom: -0.24em; }
.formula sup { top: -0.38em; }

/* 行中の数式（表のセルなど） */
.mi {
  font-family: "Cambria Math","Cambria","Times New Roman","Yu Gothic",serif;
  font-size: 11pt;
}
.mi sub, .mi sup { font-size: 68%; line-height: 0; position: relative; }
.mi sub { bottom: -0.24em; }
.mi sup { top: -0.38em; }

/* 場合分け（大きな中括弧） */
.formula.cases { display: flex; align-items: center; }
.formula.cases .fx { min-width: 0; margin-right: 6px; }
.formula.cases .brace {
  font-size: 3.1em; line-height: 0.75; font-weight: 300; margin-right: 8px;
}
.formula.cases .rows > div { line-height: 1.75; }
.formula.cases .cv { display: inline-block; min-width: 8.5em; }
.formula.cases .cc {
  font-family: "Yu Gothic","Meiryo",sans-serif; font-size: 9.5pt; color: #55636f;
}

/* ---- 表 ---- */
table {
  border-collapse: collapse; width: 100%;
  margin: 10px 0; font-size: 9pt; break-inside: avoid;
}
th, td {
  border: 1px solid #aeb9c6; padding: 4px 7px;
  text-align: left; vertical-align: top; line-height: 1.55;
}
th { background: #e8eef5; color: #17334f; font-weight: bold; }
tr:nth-child(even) td { background: #fafbfc; }

/* ---- 引用（注記ボックス） ---- */
blockquote {
  margin: 10px 0; padding: 7px 12px;
  border-left: 5px solid #4a9a78; background: #f2faf6; border-radius: 3px;
}
blockquote p { margin: 3px 0; }

ul, ol { margin: 7px 0; padding-left: 1.6em; }
li { margin: 3px 0; }

hr { border: none; border-top: 1px solid #dfe3e8; margin: 14px 0; }

a { color: #1d4266; text-decoration: none; }

/* 見出し直後の孤立行を防ぐ */
p, li { orphans: 2; widows: 2; }
"""

doc = (
    "<!DOCTYPE html>\n<html lang=\"ja\"><head><meta charset=\"utf-8\">\n"
    "<title>" + title + "</title>\n<style>" + CSS + "</style>\n</head><body>\n"
    + html_body + "\n</body></html>\n"
)

with io.open(dst, "w", encoding="utf-8") as f:
    f.write(doc)

print("html written:", dst, len(doc), "bytes")
