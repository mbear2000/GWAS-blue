GWAS-submit BLUE 扩展：本地测试版
=====================================

目的
----
这是一套“旁路测试版”，不直接覆盖你现有的 server.py 和 workflow_v2.sh。

单表型：
    仍然调用原来的 server.py / workflow_v2.sh 逻辑。

2个及以上表型：
    自动切换到 Multi-environment BLUE + GWAS。
    先直接生成 lme4 输入文件，再计算 BLUE，随后把 BLUE 结果交给
    原 workflow_v2.sh 的 EMMAx 转换、manifest、GWAS、绘图、peak SNP、
    归档流程。

文件
----
server_blue.py
    在现有 server.py 上增加 BLUE 多文件逻辑。

workflow_blue_v1.sh
    BLUE prepare 包装器。
    单表型直接 exec 原 workflow_v2.sh。
    BLUE 的 admin/final 同样 exec 原 workflow_v2.sh。

web/app.js
    测试版网页脚本。
    增加 Add phenotype、多个输入文件、每文件环境前缀、BLUE label、
    BLUE output 和预览。

start_blue.cmd
    使用 8766 端口启动本地测试版，不和正式 8765 端口混用。

安装 / 测试
-----------
1. 建议先复制一份你当前 GitHub 项目目录作为测试目录。

2. 把本压缩包中的：
       server_blue.py
       workflow_blue_v1.sh
       start_blue.cmd
   放到项目根目录，和 server.py、workflow_v2.sh 同级。

3. 用本压缩包中的：
       web/app.js
   替换测试目录里的：
       web/app.js

   正式目录建议暂时不要替换。

4. 双击：
       start_blue.cmd

5. 浏览器打开：
       http://127.0.0.1:8766/

6. 测试单表型：
       只添加一个文件。
   页面应显示：
       Single-phenotype GWAS
   此时后端会走原 workflow_v2.sh。

7. 测试 BLUE：
       添加两个或更多文件。
   页面应显示：
       Multi-environment BLUE + GWAS

   例如：
       GPallInd_LSh2022_549.csv
       GPallInd_LSh2023_549.csv
       GPallInd_LSh2024_549.csv

   自动前缀通常会识别：
       LSh2022_549_
       LSh2023_549_
       LSh2024_549_

   BLUE label：
       LSh-blue

   BLUE output：
       GPallInd_LSh_blue

服务器 BLUE 流程
----------------
上传的多个原始表型会发布到：
    <project>/phenotype/

直接生成：
    <project>/lme4/<原文件名>_lme4.csv

例如原列名：
    Accession,LSh2022_549_IRGP_days,LSh2022_549_IRGP_height

lme4 列名：
    Accession,IRGP_days,IRGP_height

然后调用：
    /public/home/yzhao/tool/R-4.1.2/bin/Rscript \
      <project>/program/generate_blue_multi_env_allTrait.R \
      <project>/lme4/GPallInd_LSh_blue.txt \
      LSh-blue \
      <多个 *_lme4.csv>

BLUE 程序来源：
    /data9/home/yzhao/program/EMMAx/generate_blue_multi_env_allTrait.R

并复制到：
    <project>/program/generate_blue_multi_env_allTrait.R

随后额外保留：
    <project>/lme4/GPallInd_LSh_blue.csv

之后把 BLUE txt 交给原 workflow_v2.sh prepare。
从 EMMAx 输入转换、manifest、GWAS、绘图、peak SNP 和归档开始，
继续复用你当前 GitHub 版本。

覆盖检查
--------
BLUE 模式启动 SecureCRT 后，会一次性检查：

    phenotype/<所有上传的原始表型>
    lme4/<BLUE output>.txt
    lme4/<BLUE output>.csv
    phenotype/<BLUE output>.txt

只要其中已有目标，就进入原网页的覆盖确认机制。
用户选择“继续并覆盖”后，才进入 prepare。

重要说明
--------
1. 这是第一版本地/服务器联调版，故意没有覆盖正式 server.py/workflow_v2.sh。
   等你实际测试 BLUE 程序的输出表头和目录结果后，再把逻辑合并成正式完整版。

2. 当前 server.py 的 HTTP 请求上限约 29 MB。
   因此 BLUE 多文件测试时，所有文件 base64 后的总请求不要太大。
   普通表型文件通常远低于这个限制。
   正式版可以把该限制改成按文件数量计算。

3. BLUE 输入文件的去前缀规则是“只删除性状名称开头完全匹配的前缀”，
   不做全局替换，因此不会误删性状名中间的年份或字符串。

4. 多个环境去前缀后：
       性状名称
       性状数量
       性状顺序
   必须完全相同，否则在本地预览/提交前停止。

5. 如果真实 generate_blue_multi_env_allTrait.R 的输出第一列、分隔符或表头
   与预期特殊不同，把一次测试日志和生成的 BLUE txt 发给我，我再按真实输出
   调整正式版。
