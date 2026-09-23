#Requires -Version 5.1
<#
  build_pdf.ps1
    技術解説書.md から 技術解説書.pdf を作る。

    通常は親フォルダの make_pdf.bat から呼ばれる。
    直接実行することもできる:
        powershell -ExecutionPolicy Bypass -File pdf_tools\build_pdf.ps1

    別のファイルを変換したい場合:
        ... -File pdf_tools\build_pdf.ps1 -Md 入力.md -Out 出力.pdf

  処理の流れ
    1. md2html.py  : Markdown -> 印刷用 HTML（数式は Cambria Math で組む）
    2. Chrome      : HTML -> PDF（--headless --print-to-pdf）
    3. stamp.py    : ページ番号とメタデータを刷り込む

  必要なもの
    ・Python 3 と markdown / pypdf / reportlab
        pip install markdown pypdf reportlab
    ・Google Chrome （見つからなければ Microsoft Edge を使う）
#>
param(
    [string]$Md  = "",
    [string]$Out = ""
)

# Chrome は起動時に無害な警告を stderr へ出す。Windows PowerShell 5.1 では
# それが NativeCommandError になってしまうので Stop にしない。
# 代わりに各工程の結果を明示的に検査する。
$ErrorActionPreference = "Continue"

$sp   = Split-Path -Parent $MyInvocation.MyCommand.Path   # pdf_tools
$proj = Split-Path -Parent $sp                            # プロジェクト直下

if (-not $Md)  { $Md  = Join-Path $proj "技術解説書.md" }
if (-not $Out) { $Out = Join-Path $proj "技術解説書.pdf" }

$html = Join-Path $sp "_work.html"      # 中間ファイル（毎回上書き）
$raw  = Join-Path $sp "_work.pdf"       # 中間ファイル（ページ番号なし）

function Fail($msg) {
    Write-Host ""
    Write-Host "[ERROR] $msg" -ForegroundColor Red
    exit 1
}

if (-not (Test-Path $Md)) { Fail "入力ファイルが見つかりません: $Md" }

# --- Python の確認 --------------------------------------------------
$py = (Get-Command python -ErrorAction SilentlyContinue)
if (-not $py) { $py = (Get-Command py -ErrorAction SilentlyContinue) }
if (-not $py) { Fail "python が見つかりません。Python 3 をインストールしてください。" }

& $py.Source -c "import markdown, pypdf, reportlab" 2>&1 | Out-Null
if ($LASTEXITCODE -ne 0) {
    Fail "必要な Python パッケージが不足しています。次を実行してください:`n        pip install markdown pypdf reportlab"
}

# --- Chrome / Edge を探す -------------------------------------------
$browser = $null
foreach ($p in @(
    "$env:ProgramFiles\Google\Chrome\Application\chrome.exe",
    "${env:ProgramFiles(x86)}\Google\Chrome\Application\chrome.exe",
    "$env:ProgramFiles\Microsoft\Edge\Application\msedge.exe",
    "${env:ProgramFiles(x86)}\Microsoft\Edge\Application\msedge.exe")) {
    if (Test-Path $p) { $browser = $p; break }
}
if (-not $browser) { Fail "Chrome も Edge も見つかりませんでした。" }

# --- 1. Markdown -> HTML --------------------------------------------
Write-Host "[1/3] Markdown -> HTML ..."
Remove-Item $html -Force -ErrorAction SilentlyContinue
& $py.Source (Join-Path $sp "md2html.py") $Md $html `
    "MPS法シミュレーションとDirectX11可視化 技術解説書" | Out-Null
if ($LASTEXITCODE -ne 0 -or -not (Test-Path $html)) { Fail "md2html.py に失敗しました。" }

# --- 2. HTML -> PDF --------------------------------------------------
Write-Host "[2/3] HTML -> PDF ($(Split-Path -Leaf $browser)) ..."
Remove-Item $raw -Force -ErrorAction SilentlyContinue
$url = "file:///" + ($html -replace '\\', '/')
& $browser --headless=new --disable-gpu --no-sandbox --no-pdf-header-footer `
           --virtual-time-budget=20000 --print-to-pdf="$raw" $url | Out-Null

# Chrome はプロセス終了後も非同期に書き込むことがある。
# サイズが 2 回続けて同じになるまで待ってから次へ進む。
$prev = -1
for ($i = 0; $i -lt 90; $i++) {
    Start-Sleep -Milliseconds 500
    if (-not (Test-Path $raw)) { continue }
    $len = (Get-Item $raw).Length
    if ($len -gt 0 -and $len -eq $prev) { break }
    $prev = $len
}
if (-not (Test-Path $raw)) { Fail "ブラウザが PDF を出力しませんでした。" }
if ((Get-Item $raw).LastWriteTime -lt (Get-Item $html).LastWriteTime) {
    Fail "PDF が HTML より古いままです（書き込みが完了していません）。"
}

# --- 3. ページ番号とメタデータ ---------------------------------------
Write-Host "[3/3] ページ番号を付与 ..."
& $py.Source (Join-Path $sp "stamp.py") $raw $Out | Out-Null
if ($LASTEXITCODE -ne 0 -or -not (Test-Path $Out)) { Fail "stamp.py に失敗しました。" }

Write-Host ""
Write-Host ("完了: {0}  ({1:N0} KB)" -f $Out, ((Get-Item $Out).Length / 1KB)) -ForegroundColor Green
exit 0
