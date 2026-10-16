# 数学物理方法 AI 助教

给 [DeepSeek Harness](https://github.com/deepseek-ai/deepseek-harness)（下称 DSH）用的一条**斜杠命令**：
在任意对话里输入 `/cyber-wuhua 你的问题`，模型就按「数学与理论物理助教」的规范来答疑、推导、
讲概念、议习题——教材优先、符号一致、推导完整、公式能正确渲染。

- 菜单显示名：**Assistant（Methods of Mathematical Physics）**
- 命令名：**`/cyber-wuhua`**（小写连字符，原因见「[常见问题](#常见问题)」）
- 附带能力：**图片、Word、PowerPoint、Excel、PDF** 题目或资料都能直接上传给它读
- 随附教材：吴崇试《数学物理方法》21 章 + 林琼桂《数学物理方法》13 章，**开箱即用**

---

## 这个文件夹里有什么

| 文件 | 作用 | 必需？ |
| --- | --- | --- |
| [`command/cyber-wuhua.md`](command/cyber-wuhua.md) | **命令本体**。整条指令的提示词都在这里，拷到 DSH 命令目录即生效 | ✅ 必需 |
| [`notation.md`](notation.md) | **符号规范**。从教材提取的符号、定义、正负号约定；助教回答前会读它 | ✅ 必需 |
| [`textbook/`](textbook/) | **随附教材原文**（35 个文件），拷进工作区即可用 | ✅ 必需 |
| [`textbook/NOTICE.md`](textbook/NOTICE.md) | 第三方内容的版权与使用范围说明 | ✅ 必读 |
| [`安装说明.md`](安装说明.md) | 从零装好的完整步骤，含插件补丁与验证方法 | ✅ 必读 |
| [`使用说明.md`](使用说明.md) | 日常怎么问、附件怎么用、公式渲染注意什么、怎么维护 | ✅ 必读 |
| [`教材准备.md`](教材准备.md) | 教材命名约定；用别的教材时怎么替换 | 换教材才读 |
| [`patch/附件支持补丁.md`](patch/附件支持补丁.md) | 让命令能收附件的插件补丁：确切改哪两处、怎么改、怎么回退 | 需要附件才装 |
| [`LICENSE`](LICENSE) | MIT；**仅覆盖本项目原创部分**，不含 `textbook/` | — |
| [`verify/检查安装.sh`](verify/检查安装.sh) · [`verify/检查安装.ps1`](verify/检查安装.ps1) | 一条命令检查装得对不对 | 建议 |

> ⚠️ `textbook/` 是第三方教材内容，版权归原作者与出版方，**不适用根目录的 MIT 许可**。
> 详见 [`textbook/NOTICE.md`](textbook/NOTICE.md)。

---

## 装它需要三样东西

```
DSH 客户端
   ├── 插件：dsh-slash-command-filesystem     ← 把 .md 变成 /命令，并监听目录变化
   │      └── 一条补丁：让它能接收附件          ← patch/附件支持补丁.md
   ├── 命令本体：cyber-wuhua.md                ← command/
   └── 工作区素材：textbook/ + notation.md      ← 本包提供，拷进工作区即可
```

四样都要就位，但**只有插件那一项需要动别人的代码**。安装说明里两条路都写了：
打补丁（有附件）或原样使用（无附件）。

---

## 五分钟上手

```powershell
# 1) 装插件
dsh plugin --profile desktop add dsh-slash-command-filesystem

# 2) 放命令本体（用户级目录 → 所有工作区、所有会话都可见）
New-Item -ItemType Directory -Force "$env:USERPROFILE\.dsh\commands" | Out-Null
Copy-Item ".\command\cyber-wuhua.md" "$env:USERPROFILE\.dsh\commands\cyber-wuhua.md"

# 3) 把随附教材与符号规范铺进你的工作区
#    假设工作区是 D:\study
$ws = "D:\study"
New-Item -ItemType Directory -Force $ws | Out-Null
Copy-Item ".\textbook\*.md"  $ws -Force      # 34 章 + contents.md
Copy-Item ".\notation.md"    $ws -Force      # 符号规范

# 4) 检查（Windows 默认禁止运行 .ps1，需带 ExecutionPolicy Bypass）
powershell -NoProfile -ExecutionPolicy Bypass -File ".\verify\检查安装.ps1"
```

然后在 DSH 里**把工作区切到 `D:\study`**，输入框敲 `/cyber` —— 菜单里出现即可用。

> `--profile desktop` 是桌面客户端的 profile 名；CLI／网页端请换成你自己的 profile 名
> （`dsh plugin list` 可以查）。
>
> 脚本编码：`检查安装.ps1` 带 UTF-8 BOM（Windows PowerShell 5.1 才能正确读中文），
> `检查安装.sh` 是 LF + UTF-8（无 BOM）。**不要把两者的编码对调**。

---

## 设计取舍（开源时想让你知道的）

**教材与提示词是解耦的。** 命令的提示词只依赖**文件命名约定**（`chapter-NN-书名.md`）与
**你自己的 `notation.md`**，不含任何教材正文。所以换一本教材只需改提示词里两行文件名模式、
重写 `notation.md`——其余部分与具体教材无关。见 [`教材准备.md`](教材准备.md)。

**`notation.md` 是「可替换件」，也是整个设计最关键的一处解耦。** 本包那份是从吴崇试、林琼桂
两套教材**逐条核对原文**得出的，助教把它当作本项目的默认符号规范，优先级高于模型自己的习惯。
换教材就重写它，助教行为会自动跟着变。

**提示词写的是「通用助教规范」而不是「某本书的答案」**：教材优先、不得擅自替换教材符号、
新符号先定义、两种约定冲突时必须声明、以帮助理解为目标、超出教材内容要标注、
偏导保持哪些变量不变、区分有限变化/微分/变分、区分充分必要充要、明确近似阶数、
注意定义域与适用范围……这些对任何一本数学物理方法教材都成立。

**KaTeX 那一段是踩坑换来的。** 提示词里专门有一节讲 `$$` 必须独占一行、
不能两个独立公式首尾相接、不能用 `$$$$`。这不是洁癖：DSH 前端按 `$$` 切分数学块，
粘连会切出错位的块，KaTeX 就报 `Expected 'EOF, got ...`，整段公式变红。
详见[使用说明](使用说明.md#公式渲染)。

**附件补丁改的是别人的文件。** `dsh-slash-command-filesystem` 有两处没实现对附件的支持
（注册时没声明 `input.attachments`；处理器没收 `invocation.attachments`）。
补丁就两处、约 10 行，但**它打在 `node_modules` 里，重新安装插件会被覆盖**——
补丁说明里写了怎么重打和怎么回退。

---

## 常见问题

**为什么不是 `/Cyber WuHua`？**
DSH 命令名文法只允许 `^[a-z][a-z0-9_-]*$`，**不接受大写和空格**；而且解析器只把第一个空格前的内容
当命令名，`/Cyber WuHua 问题` 会被当成「命令 `Cyber` + 参数」，查不到就整行退化成普通消息。
`/cyber-wuhua` 是唯一既能被解析、又能被菜单命中的写法。

**为什么命令「找不到」了？**
先说结论：**别 fork 这个插件**。把插件复制到别处再 link 进去，插件会因无法解析
`@deepseek-ai/*` 这类由宿主提供的裸导入而**加载失败**，命令随之消失。
正确做法是就地给已安装的插件打补丁——见 [`patch/附件支持补丁.md`](patch/附件支持补丁.md)。

**上传图片报 `Declared image type does not match its bytes`？**
文件扩展名和真实字节不符（常见于从聊天软件另存、浏览器「另存为」类型选错）。
用 Pillow 按真实格式重新保存，或改扩展名。可参考[使用说明](使用说明.md#附件)里的排查方法。

**公式变红 `KaTeX parse error`？**
见[使用说明 · 公式渲染](使用说明.md#公式渲染)。

**能把教材换成我自己的吗？**
可以，而且推荐。见 [`教材准备.md`](教材准备.md)。

---

## 许可与出处

| 部分 | 许可 / 版权 |
| --- | --- |
| 命令提示词、`notation.md`、文档、校验脚本 | **MIT**（见 [`LICENSE`](LICENSE)） |
| `textbook/` 教材原文 | **不属于本项目**，版权归原作者与出版方，见 [`textbook/NOTICE.md`](textbook/NOTICE.md) |
| 插件 `dsh-slash-command-filesystem` | MIT，作者 [rengy](https://github.com/rengy/dsh-slash-command-filesystem)；本包只提供补丁说明，不分发其代码 |
| DSH 本体 | DeepSeek 官方项目 |

> 若要再分发本包，请一并保留 `LICENSE` 与 `textbook/NOTICE.md`。
> 若无法确认教材内容的再分发合规性，**删掉 `textbook/` 目录再分发即可**——
> 助教不依赖特定教材。
