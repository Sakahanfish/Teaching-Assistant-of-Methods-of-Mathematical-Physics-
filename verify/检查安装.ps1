#!/usr/bin/env pwsh
# 检查「数学物理方法 AI 助教」安装是否正确。
# 用法：  .\检查安装.ps1
#         .\检查安装.ps1 -Profile web

param(
    [string]$Profile = "desktop",
    [string]$DshHome = $env:DSH_HOME
)

$ErrorActionPreference = 'Continue'
if ([string]::IsNullOrWhiteSpace($DshHome)) {
    $DshHome = Join-Path $env:USERPROFILE '.dsh'
}

$pass = 0
$fail = 0
$warn = 0

function Ok($msg)   { Write-Host "  [OK]   $msg" -ForegroundColor Green;  $script:pass++ }
function Bad($msg)  { Write-Host "  [FAIL] $msg" -ForegroundColor Red;    $script:fail++ }
function Warn($msg) { Write-Host "  [WARN] $msg" -ForegroundColor Yellow; $script:warn++ }
function Head($msg) { Write-Host "`n== $msg" -ForegroundColor Cyan }

Write-Host "DSH_HOME   : $DshHome"
Write-Host "profile    : $Profile"

# ---------------------------------------------------------------- 1. 命令本体
Head "1. 命令本体"

$cmd = Join-Path $DshHome 'commands\cyber-wuhua.md'
if (-not (Test-Path $cmd)) {
    Bad "找不到 $cmd —— 请把 command/cyber-wuhua.md 拷到该位置"
} else {
    $raw = [System.IO.File]::ReadAllText($cmd)
    Ok "命令文件存在（$($raw.Length) 字符）"

    if ($raw -match '(?s)^---\r?\n(.*?)\r?\n---') {
        Ok "frontmatter 块存在"
        $fm = $Matches[1]
        foreach ($key in @('description', 'argument-hint')) {
            if ($fm -match "(?m)^$key\s*:\s*\S") { Ok "frontmatter 含 $key" }
            else { Bad "frontmatter 缺少 $key" }
        }
    } else {
        Bad "frontmatter 块缺失或格式错误（该文档会被跳过）"
    }

    if ($raw -match '\$ARGUMENTS') { Ok '正文含 $ARGUMENTS 占位符' }
    else { Bad '正文缺少 $ARGUMENTS —— 用户的问题不会被传入' }

    $body = $raw.Substring($raw.IndexOf('---', 3) + 3)
    if ($body -match 'notation\.md') { Ok '正文引用了 notation.md' } else { Warn '正文未引用 notation.md' }
    if ($body -match 'KaTeX')        { Ok '正文含 KaTeX 渲染约束' } else { Warn '正文缺 KaTeX 约束（公式可能变红）' }
}

# ---------------------------------------------------------------- 2. 插件
Head "2. 插件 dsh-slash-command-filesystem"

$pluginDir = Join-Path $DshHome "profiles\$Profile\node_modules\dsh-slash-command-filesystem"
$pluginLib = Join-Path $pluginDir 'lib\index.js'
if (Test-Path $pluginLib) {
    Ok "插件已安装：$pluginDir"
} else {
    Bad "插件未安装。请运行： dsh plugin --profile $Profile add dsh-slash-command-filesystem"
}

$profilePkg = Join-Path $DshHome "profiles\$Profile\package.json"
if (Test-Path $profilePkg) {
    try {
        $pkg = Get-Content $profilePkg -Raw | ConvertFrom-Json
        $bundles = @($pkg.dsh.profile.bundles)
        if ($bundles -contains 'dsh-slash-command-filesystem') {
            Ok 'profile bundle 列表含该插件'
        } else {
            Bad "profile bundle 列表不含该插件（当前：$($bundles -join ', ')）"
        }
    } catch {
        Warn "无法解析 $profilePkg"
    }
} else {
    Warn "找不到 $profilePkg"
}

# ---------------------------------------------------------------- 3. 附件补丁
Head "3. 附件支持补丁（可选）"

if (Test-Path $pluginLib) {
    $src = [System.IO.File]::ReadAllText($pluginLib)
    $hasGate = $src -match 'attachments:\s*true'
    $hasForward = $src -match 'invocation\.attachments'

    if ($hasGate)    { Ok '准入关：已声明 input.attachments' }
    else             { Warn '准入关未通过 —— 附件提交会被拒（见 patch/附件支持补丁.md）' }
    if ($hasForward) { Ok '投递关：已转发 invocation.attachments' }
    else             { Warn '投递关未通过 —— 附件会被静默丢弃（见 patch/附件支持补丁.md）' }

    if ($hasGate -and $hasForward) {
        Write-Host "         注意：补丁在启动时加载，改完需重启一次 DSH 才生效。" -ForegroundColor DarkGray
    }
} else {
    Warn '插件未安装，跳过补丁检查'
}

# ---------------------------------------------------------------- 4. 工作区素材
Head "4. 工作区素材（在命令被调用的工作区里检查）"

$cwd = (Get-Location).Path
Write-Host "          当前目录：$cwd" -ForegroundColor DarkGray

if (Test-Path (Join-Path $cwd 'notation.md')) {
    Ok 'notation.md 在当前工作区（助教会读它作为符号规范）'
} else {
    Warn '当前工作区没有 notation.md —— 助教会跳过这一步并说明；建议放一份'
}

$chapters = @(Get-ChildItem -Path $cwd -Filter 'chapter-*-*.md' -File -ErrorAction SilentlyContinue)
if ($chapters.Count -gt 0) {
    $books = $chapters | ForEach-Object {
        if ($_.BaseName -match '^chapter-\d+-(.+)$') { $Matches[1] }
    } | Sort-Object -Unique
    Ok "找到 $($chapters.Count) 个教材章节文件，教材标识：$($books -join ', ')"
} else {
    Warn '当前工作区没找到 chapter-*-*.md 教材文件 —— 助教将无法引用原文（见 教材准备.md）'
}

# ---------------------------------------------------------------- 汇总
Head "汇总"
Write-Host "  OK: $pass    WARN: $warn    FAIL: $fail"
if ($fail -eq 0) {
    Write-Host "`n  可以用了：在输入框敲 /cyber 看看菜单里有没有它。" -ForegroundColor Green
    exit 0
} else {
    Write-Host "`n  有 $fail 项失败，请按上面的提示处理；详见 安装说明.md。" -ForegroundColor Red
    exit 1
}
