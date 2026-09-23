# -*- coding: utf-8 -*-
"""PDF にページ番号（N / total）を刷り込み、メタデータを設定する"""
import sys, io
from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen import canvas

src, dst = sys.argv[1], sys.argv[2]

reader = PdfReader(src)
total = len(reader.pages)
writer = PdfWriter()

for i, page in enumerate(reader.pages, start=1):
    w = float(page.mediabox.width)
    h = float(page.mediabox.height)

    # 1 ページ分のオーバーレイをメモリ上に作る
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=(w, h))
    c.setFont("Helvetica", 8.5)
    c.setFillGray(0.45)
    # 数字だけなので欧文フォントで足りる（日本語フォント埋め込み不要）
    c.drawCentredString(w / 2.0, 26.0, "%d / %d" % (i, total))
    c.save()
    buf.seek(0)

    overlay = PdfReader(buf).pages[0]
    page.merge_page(overlay)
    writer.add_page(page)

writer.add_metadata({
    "/Title":    "MPS法シミュレーションとDirectX 11可視化 技術解説書",
    "/Subject":  "粒子法(MPS法)による流体シミュレーションとDirectX 11による可視化の技術解説",
    "/Creator":  "markdown -> HTML -> Chrome print-to-pdf",
})

with open(dst, "wb") as f:
    writer.write(f)

print("stamped %d pages -> %s" % (total, dst))
