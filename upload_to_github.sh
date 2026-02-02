#!/bin/bash
# =============================================================================
# ProteomicsApp GitHub上传脚本
# =============================================================================

echo "============================================================================="
echo "ProteomicsApp GitHub 上传助手"
echo "============================================================================="
echo ""

# 检查是否在ProteomicsApp目录中
if [ ! -f "DESCRIPTION" ]; then
    echo "错误: 请在ProteomicsApp目录中运行此脚本"
    exit 1
fi

# 检查git是否已初始化
if [ ! -d .git ]; then
    echo "初始化Git仓库..."
    git init
    echo "✓ Git仓库已初始化"
    echo ""
fi

# 添加所有文件
echo "步骤 1/4: 添加文件到Git..."
git add .
echo "✓ 文件已添加"
echo ""

# 显示状态
echo "当前状态:"
git status --short
echo ""

# 提交
echo "步骤 2/4: 提交更改..."
read -p "请输入提交信息 (默认: 'Initial commit'): " commit_msg
commit_msg=${commit_msg:-"Initial commit"}
git commit -m "$commit_msg"
echo "✓ 更改已提交"
echo ""

# 询问GitHub仓库URL
echo "步骤 3/4: 设置远程仓库..."
echo ""
echo "请先在GitHub上创建新仓库（如果尚未创建）："
echo "  1. 访问 https://github.com/new"
echo "  2. 仓库名称: ProteomicsApp"
echo "  3. 描述: A comprehensive proteomics analysis platform"
echo "  4. 选择 Public 或 Private"
echo "  5. 不要勾选 'Initialize with README'（我们已经有了）"
echo "  6. 点击 'Create repository'"
echo ""
read -p "请输入GitHub仓库URL (如: https://github.com/username/ProteomicsApp.git): " repo_url

if [ -z "$repo_url" ]; then
    echo "错误: 仓库URL不能为空"
    exit 1
fi

# 检查remote是否已存在
if git remote | grep -q "origin"; then
    echo "更新远程仓库URL..."
    git remote set-url origin "$repo_url"
else
    echo "添加远程仓库..."
    git remote add origin "$repo_url"
fi
echo "✓ 远程仓库已设置"
echo ""

# 推送
echo "步骤 4/4: 推送到GitHub..."
git branch -M main
git push -u origin main
echo "✓ 代码已推送到GitHub"
echo ""

# 完成
echo "============================================================================="
echo "上传完成！"
echo "============================================================================="
echo ""
echo "您的仓库地址: $repo_url"
echo ""
echo "下一步操作："
echo "  1. 访问仓库页面，检查文件是否正确上传"
echo "  2. 更新 README.md 中的 GitHub 链接"
echo "  3. （可选）创建 Release 标签"
echo ""
echo "团队成员现在可以使用以下命令安装："
echo "  devtools::install_github(\"your-username/ProteomicsApp\")"
echo ""
