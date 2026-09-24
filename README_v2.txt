GWAS-submit BLUE 本地测试版 v2
================================

本版修复 v1 的关键问题
---------------------
v1 在远端把 BLUE 包装脚本上传成：
    .gwas-web/<runid>/workflow.sh

但包装脚本随后需要：
    .gwas-web/<runid>/workflow_v2.sh

v1 没有把原始 workflow_v2.sh 一起上传，因此 prepare 会在建立目标项目目录之前退出。
因为退出发生在原 workflow_v2.sh 的日志/exit-marker 初始化之前，网页会一直显示 prepare。

v2 修复为：
1. SecureCRT 上传 BLUE wrapper 之后，再显式上传原 workflow_v2.sh。
2. BLUE wrapper 从启动第一行就追加写入：
       .gwas-web/<runid>/run.log
3. BLUE wrapper 自己发生错误时也会写：
       .gwas-web/<runid>/prepare.exit
   因此网页能立即显示失败，而不会无限等待。
4. start_blue.cmd 会自动打开：
       http://127.0.0.1:8766/

安装
----
在你当前 GWAS-submit 的“测试副本”中保证已有：
    server.py
    workflow_v2.sh
    web/index.html

把 v2 包中的以下文件复制进去：
    server_blue.py
    workflow_blue_v2.sh
    start_blue.cmd

并用：
    web/app.js
替换测试副本中的：
    web/app.js

不要删除原 workflow_v2.sh。

然后双击：
    start_blue.cmd

旧的 v1 卡住任务
----------------
旧任务没有进入原 workflow_v2.sh prepare，也没有提交 qsub GWAS 作业。
从截图表现看，它停在 BLUE wrapper 找不到 workflow_v2.sh 的位置。

可以关闭旧 SecureCRT 测试标签页，不要对旧任务点击“恢复”。
使用 v2 重新发起一个新的测试任务。

服务器快速验证
--------------
v2 启动 BLUE 后，fat2 上应该很快能看到：

    /data9/home/yzhao/GWAS_IRGSP1.0/.gwas-web/<runid>/workflow.sh
    /data9/home/yzhao/GWAS_IRGSP1.0/.gwas-web/<runid>/workflow_v2.sh
    /data9/home/yzhao/GWAS_IRGSP1.0/.gwas-web/<runid>/run.log

随后目标项目目录应该建立，例如：

    /data9/home/yzhao/GWAS_IRGSP1.0/GPallJapmiss20_UAV_549lines_5years_Re-fittingBlue/

至少会很快出现：
    phenotype/
    lme4/

如果 BLUE R 程序或输入有问题，网页现在会明确进入失败/需要处理状态，
run.log 中会保留实际 ERROR，不会再永久停在 prepare。

建议第一次重新测试时，在服务器另开终端观察：

    ls -lt /data9/home/yzhao/GWAS_IRGSP1.0/.gwas-web/ | head

找到最新 runid 后：

    tail -f /data9/home/yzhao/GWAS_IRGSP1.0/.gwas-web/<runid>/run.log

这样 BLUE 的 lme4 转换和 Rscript 输出都会实时看到。
