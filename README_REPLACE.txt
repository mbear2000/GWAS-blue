GWAS-submit BLUE 完整替换版 v9.2
==============================

本版主要统一 BLUE 命名，不改变已经跑通的 BLUE / GWAS 核心分析逻辑。

例如：
    population = GPallInd
    本次标签 = LSH-blue

文件/目录名中的短横线转换为下划线，统一使用：
    GPallInd_LSH_blue

远程结果命名：
    phenotype/GPallInd_LSH_blue.txt
    EMMAx.Result/hIBS/GPallInd_LSH_blue/
    sh.file/emmax/GPallInd_LSH_blue/
    sh.file/plotGWAS/GPallInd_LSH_blue/
    sigSNP/GPallInd_LSH_blue.list
    sigSNP/GPallInd_LSH_blue_geneList.csv

不再出现：
    GPallInd_GPallInd_LSH_blue_BLUE_multi.list
    GPallInd_LSH_blue_BLUE_multi_20260923-1913

GWAS trait 仍保留用户填写标签中的短横线，例如：
    GPallInd_LSH-blue_IRGP_days

网页运行历史显示：
    GPallInd_LSH_blue

下载日志文件名：
    GPallInd_LSH_blue.log

浏览器允许通知后：
    分析完成 -> Windows 通知
    失败/需要处理 -> Windows 通知

workflow_v2.sh 和 workflow_blue_v8.sh 必须同时存在。
不要用 workflow_blue_v8.sh 覆盖 workflow_v2.sh。

启动：
    start_blue.cmd

新任务日志开头应看到：
    BLUE_WRAPPER_V8_START


七、关闭本地 BLUE 服务
---------------------

新增：

    stop_blue.cmd

以后不需要在 start_blue.cmd 窗口里按 Ctrl+C。

双击 stop_blue.cmd 后，它只检查：

    127.0.0.1:8766

并且只有确认该端口对应的进程命令行包含：

    server_blue.py

才会结束该进程。

它不会：

    - 停止 8765 的原 GWAS 服务
    - 删除本地 runs 历史
    - 删除远程 .gwas-active
    - 取消 qstat 中的远程作业
    - 修改任何 GWAS / BLUE 结果

如果 8766 已经停止，会直接提示：

    No process is listening on port 8766. GWAS BLUE is already stopped.


八、v8.2 同名覆盖确认修复
------------------------

修复：远程检查发现同名输入/输出后，SecureCRT 打印 yes 但随后通信超时，
网页进入“需要处理”，没有显示覆盖确认按钮。

原因：v8.1 的检查命令命中同名文件后执行 exit 0，导致 ExecScalar 还没来得及
输出结束标记，SecureCRT 就失去当前 shell，因此被误判为超时。

v8.2：
- 同名检查不再 exit；
- 正常返回 yes/no；
- 状态进入“等待覆盖确认”；
- 网页显示“继续并覆盖”和“停止本次分析”；
- 确认覆盖后，只覆盖本次同名输入/输出；
- 不同标签的既有结果不会删除。

同一工作目录中，以前跑过 LSH-blue，现在上传相同环境文件但标签改为
LocYear-blue 时，输入 phenotype 文件名可能仍相同，因此仍会提示一次覆盖。
这是正常的，确认继续即可。


九、v8.3 qsub 提交间隔
---------------------

每次 qsub 成功提交并记录 Job ID 后，等待时间改为：

    sleep 1

即每个 qsub 之间等待 1 秒，不再等待 5 秒。

注意：这里只修改“每次 qsub 后的短暂停顿”。
队列状态轮询逻辑仍保持原来的 120 秒 / 10 秒检查规则，不作改动。


十、v8.4 下载日志文件名
----------------------

窗口显示名保持：

    GPallInd_LocYear_blue

下载日志时保留日期时间，文件名改为：

    GPallInd_LocYear_blue_20260923-203112.log

不再使用：

    GPallInd_LocYear_blue_GWASrun_GPallInd_20260923-203112.log

日期时间直接取本次本地日志文件末尾的 YYYYMMDD-HHMMSS。


十一、v8.5 BLUE R 程序来源
-------------------------

以下两个 BLUE R 程序不再从 Windows 本地上传：

    generate_blue_multi_env_allTrait_fixed.R
    generate_blue_multiEnvYear_allTrait_fixed.R

每次 prepare 时统一从服务器复制：

    /data9/home/yzhao/program/EMMAx/generate_blue_multi_env_allTrait_fixed.R
    /data9/home/yzhao/program/EMMAx/generate_blue_multiEnvYear_allTrait_fixed.R

复制到当前项目：

    <project>/program/generate_blue_multi_env_allTrait_fixed.R
    <project>/program/generate_blue_multiEnvYear_allTrait_fixed.R

实际 Rscript 也直接调用这两个 fixed 文件。

当前仍由 Windows/本地网页上传到本次远程 staging 的程序性文件只有：

    1. workflow_blue_v8.sh
       BLUE 外层控制脚本，由 server_blue.py 指定为当前 workflow。

    2. workflow_v2.sh
       本地 workflow_v2.sh 经 BLUE 兼容补丁后上传，
       用于后续原 GWAS 主流程。

其他 GWAS Perl/R 分析程序不从 Windows 上传。
它们仍由服务器端 workflow_v2.sh 从：

    /data9/home/yzhao/GWAS_IRGSP1.0/000data_prepare/program

复制到项目/批次目录。


十二、v8.6 下载日志文件名
------------------------

窗口显示名保持原样，例如：

    GPallInd_LocYear_blue

下载日志文件名改为：

    GPallInd_LocYear_blue_20260923-203112.log

即：

    <窗口名称>_<YYYYMMDD-HHMMSS>.log

不再包含：

    _GWASrun_GPallInd


十三、v9 genomic inflation factor
---------------------------------

BLUE 和原始单表型 GWAS 都在 final 阶段最后增加 Lambda GC 计算。

程序固定来自服务器：
    /data9/home/yzhao/GWAS_IRGSP1.0/000data_prepare/program/get_genomicInflationFactor.py

不会从 Windows 本地上传。prepare 阶段现有逻辑会把 canonical program 目录
复制到本批次：
    <project>/.gwas-runs/<runid>/program/

final 阶段使用：
    <batch>/program/get_genomicInflationFactor.py

输入：
    已归档结果目录 output/ 中每个 trait 的 chr01 ps 文件：
    <result-archive>/output/<trait>_chr01_hIBS.ps

脚本会自动扩展 chr01 到 chr12。

输出：
    <result-archive>/inflationFactor/<trait>.inflationFactor

Conda：
    final 阶段先 source conda.sh，再执行：
    conda activate "${GWAS_INFLATION_CONDA_ENV:-base}"

默认环境为 base。
如果 numpy/scipy 安装在其他环境，可设置：
    export GWAS_INFLATION_CONDA_ENV=<环境名>

日志：
    INFLATION_FACTOR_START
    INFLATION_FACTOR_INPUT
    INFLATION_FACTOR_OUTPUT
    INFLATION_FACTOR_DONE
    INFLATION_FACTOR_ALL_DONE

任一 trait 缺少 chr01 ps、Python 失败或输出为空，final 会失败。


十四、v9.1 进度界面与 conda
---------------------------

运行进度新增第 5 步：

    1 准备数据与提取表型
    2 GWAS 队列计算
    3 绘图与日志检查
    4 提取 peak SNP 与已知基因
    5 计算基因组膨胀系数
    ✓ 分析结束，保存运行日志

当 final 日志出现：

    ===== GENOMIC INFLATION FACTOR =====
或：
    INFLATION_FACTOR_START

网页会把第 5 步标为正在进行。

Conda 改为直接：

    conda activate

不再指定 base 或其他环境名。

为了保证非交互 shell 中 conda activate 可用，流程仍会先找到并 source conda.sh，
然后执行不带环境名的：

    conda activate

同时网页说明中的“每次 qsub 后等待5秒”改为“等待1秒”。


十五、v9.2 BLUE R 程序名修正
---------------------------

BLUE R 程序继续从服务器：

    /data9/home/yzhao/program/EMMAx/

复制，但文件名恢复为实际使用的：

    generate_blue_multi_env_allTrait.R
    generate_blue_multiEnvYear_allTrait.R

不再使用：

    generate_blue_multi_env_allTrait_fixed.R
    generate_blue_multiEnvYear_allTrait_fixed.R

Windows 本地不上传这两个 R 程序。

本 zip 不再包含外层版本文件夹。
解压后根目录直接得到：

    server_blue.py
    workflow_v2.sh
    workflow_blue_v8.sh
    start_blue.cmd
    stop_blue.cmd
    web/
    README_REPLACE.txt

可以直接解压覆盖测试目录。
