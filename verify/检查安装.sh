#!/usr/bin/env bash
# 检查「数学物理方法 AI 助教」安装是否正确。
# 用法：  bash ./检查安装.sh
#         DSH_HOME=/path/to/.dsh PROFILE=web bash ./检查安装.sh
#
# 支持 macOS / Linux / Git Bash / WSL。

set -u

DSH_HOME="${DSH_HOME:-$HOME/.dsh}"
PROFILE="${PROFILE:-desktop}"

pass=0; fail=0; warn=0
if [ -t 1 ]; then
  G=$'\033[32m'; R=$'\033[31m'; Y=$'\033[33m'; C=$'\033[36m'; D=$'\033[90m'; N=$'\033[0m'
else
  G=''; R=''; Y=''; C=''; D=''; N=''
fi
ok()   { printf '  %s[OK]%s   %s\n'   "$G" "$N" "$1"; pass=$((pass+1)); }
bad()  { printf '  %s[FAIL]%s %s\n'   "$R" "$N" "$1"; fail=$((fail+1)); }
warn() { printf '  %s[WARN]%s %s\n'   "$Y" "$N" "$1"; warn=$((warn+1)); }
head() { printf '\n%s== %s%s\n' "$C" "$1" "$N"; }

echo "DSH_HOME   : $DSH_HOME"
echo "profile    : $PROFILE"

# ---------------------------------------------------------------- 1. 命令本体
head "1. 命令本体"

CMD="$DSH_HOME/commands/cyber-wuhua.md"
if [ ! -f "$CMD" ]; then
  bad "找不到 $CMD —— 请把 command/cyber-wuhua.md 拷到该位置"
else
  size=$(wc -c < "$CMD" | tr -d ' ')
  ok "命令文件存在（$size 字节）"

  if head -1 "$CMD" | grep -q '^---'; then
    ok "frontmatter 起始标记存在"
    fm=$(awk 'NR==1{next} /^---[[:space:]]*$/{exit} {print}' "$CMD")
    for key in description argument-hint; do
      if printf '%s\n' "$fm" | grep -qE "^$key[[:space:]]*:[[:space:]]*[^[:space:]]"; then
        ok "frontmatter 含 $key"
      else
        bad "frontmatter 缺少 $key"
      fi
    done
  else
    bad "frontmatter 缺失（该文档会被插件跳过）"
  fi

  if grep -q '\$ARGUMENTS' "$CMD"; then
    ok '正文含 $ARGUMENTS 占位符'
  else
    bad '正文缺少 $ARGUMENTS —— 用户的问题不会被传入'
  fi

  grep -q 'notation\.md' "$CMD" && ok '正文引用了 notation.md' || warn '正文未引用 notation.md'
  grep -q 'KaTeX'        "$CMD" && ok '正文含 KaTeX 渲染约束' || warn '正文缺 KaTeX 约束（公式可能变红）'
fi

# ---------------------------------------------------------------- 2. 插件
head "2. 插件 dsh-slash-command-filesystem"

PLUGIN_DIR="$DSH_HOME/profiles/$PROFILE/node_modules/dsh-slash-command-filesystem"
PLUGIN_LIB="$PLUGIN_DIR/lib/index.js"
if [ -f "$PLUGIN_LIB" ]; then
  ok "插件已安装：$PLUGIN_DIR"
else
  bad "插件未安装。请运行： dsh plugin --profile $PROFILE add dsh-slash-command-filesystem"
fi

PROFILE_PKG="$DSH_HOME/profiles/$PROFILE/package.json"
if [ -f "$PROFILE_PKG" ]; then
  if grep -q 'dsh-slash-command-filesystem' "$PROFILE_PKG"; then
    ok 'profile bundle 列表含该插件'
  else
    bad 'profile bundle 列表不含该插件'
  fi
else
  warn "找不到 $PROFILE_PKG"
fi

# ---------------------------------------------------------------- 3. 附件补丁
head "3. 附件支持补丁（可选）"

if [ -f "$PLUGIN_LIB" ]; then
  if grep -qE 'attachments:[[:space:]]*true' "$PLUGIN_LIB"; then
    ok '准入关：已声明 input.attachments'
  else
    warn '准入关未通过 —— 附件提交会被拒（见 patch/附件支持补丁.md）'
  fi
  if grep -q 'invocation\.attachments' "$PLUGIN_LIB"; then
    ok '投递关：已转发 invocation.attachments'
  else
    warn '投递关未通过 —— 附件会被静默丢弃（见 patch/附件支持补丁.md）'
  fi
  printf '  %s注意：补丁在启动时加载，改完需重启一次 DSH 才生效。%s\n' "$D" "$N"
else
  warn '插件未安装，跳过补丁检查'
fi

# ---------------------------------------------------------------- 4. 工作区素材
head "4. 工作区素材（在命令被调用的工作区里检查）"

CWD="$(pwd)"
printf '  %s当前目录：%s%s\n' "$D" "$CWD" "$N"

if [ -f "$CWD/notation.md" ]; then
  ok 'notation.md 在当前工作区（助教会读它作为符号规范）'
else
  warn '当前工作区没有 notation.md —— 助教会跳过这一步并说明；建议放一份'
fi

chapters=$(find "$CWD" -maxdepth 1 -name 'chapter-*-*.md' -type f 2>/dev/null | wc -l | tr -d ' ')
if [ "$chapters" -gt 0 ]; then
  books=$(find "$CWD" -maxdepth 1 -name 'chapter-*-*.md' -type f -exec basename {} .md \; 2>/dev/null \
          | sed -E 's/^chapter-[0-9]+-//' | sort -u | tr '\n' ' ')
  ok "找到 $chapters 个教材章节文件，教材标识：$books"
else
  warn '当前工作区没找到 chapter-*-*.md 教材文件 —— 助教将无法引用原文（见 教材准备.md）'
fi

# ---------------------------------------------------------------- 汇总
head "汇总"
printf '  OK: %d    WARN: %d    FAIL: %d\n' "$pass" "$warn" "$fail"
if [ "$fail" -eq 0 ]; then
  printf '\n  %s可以用了：在输入框敲 /cyber 看看菜单里有没有它。%s\n' "$G" "$N"
  exit 0
else
  printf '\n  %s有 %d 项失败，请按上面的提示处理；详见 安装说明.md。%s\n' "$R" "$fail" "$N"
  exit 1
fi
