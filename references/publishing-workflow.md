# 分享 Skill 到 GitHub 工作流

## 前置条件

- GitHub 账号
- 本地安装 git + gh (GitHub CLI)

## 安装 gh

```bash
# 推荐：从 GitHub Releases 直接下载（绕过 brew 网络问题）
# 1. 查最新版本号
curl -fsSL "https://api.github.com/repos/cli/cli/releases/latest" | \
  grep -oE '"browser_download_url": "[^"]*macOS[^"]*arm64[^"]*"'

# 2. 下载对应版本的 macOS arm64 zip
curl -fsSL -o /tmp/gh.zip \
  "https://github.com/cli/cli/releases/download/v2.96.0/gh_2.96.0_macOS_arm64.zip"

# 3. 解压并安装
cd /tmp && unzip -o gh.zip
sudo cp gh_2.96.0_macOS_arm64/bin/gh /usr/local/bin/
```

## 登录 GitHub

### 方式一：设备码登录（推荐）

```bash
gh auth login --web -h github.com
# 浏览器打开 https://github.com/login/device
# 输入终端显示的设备代码
```

### 方式二：钥匙串已有 Token 直接登录（备用）

如果 macOS 钥匙串中已经有 GitHub Personal Access Token，可以直接取用：

```bash
# 从钥匙串提取 Token
GH_TOKEN=$(security find-internet-password -s github.com -w 2>/dev/null)

# 用 Token 登录 gh
echo "$GH_TOKEN" | gh auth login --with-token

# 验证
gh auth status
```

⚠️ **注意**：设备码方式（--web）在 Hermes 的 PTY/后台模式下可能因回调超时而失败。如果试了多次不成功，直接用钥匙串 Token 是更稳妥的方案。

## 创建仓库并推送 Skill

```bash
# 1. 在临时目录初始化 git 仓库
mkdir -p /tmp/skill-publish
cp -r ~/.hermes/skills/research/drug-target-mr/* /tmp/skill-publish/
cd /tmp/skill-publish

# 2. 创建 README
cat > README.md << 'EOF'
# drug-target-mr — 药靶孟德尔随机化 Skill

Drug-target Mendelian Randomization 完整工作流。
涵盖人群匹配、GTEx eQTL 处理、FinnGen 结局提取、MR 分析、基因拆分等。
EOF

# 3. 初始化并推送
git init
git add -A
git commit -m "Initial commit: drug-target-mr skill"
gh repo create drug-target-mr-skill --public --push --source=.
```

## 安装已发布的 Skill

其他人可以用以下命令安装：

```bash
# 方式一：直接 URL 安装
hermes skills install \
  https://raw.githubusercontent.com/<用户名>/drug-target-mr-skill/main/SKILL.md

# 方式二：添加为 tap 源
hermes skills tap add https://github.com/<用户名>/drug-target-mr-skill
```
