GWAS-submit BLUE 完整替换版 v9.9
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


十六、v9.3 conda 激活修复
------------------------

不再查找固定 conda.sh 路径。

膨胀系数最后一步通过交互式 bash 启动：

    bash -ic

在该 shell 中直接执行：

    conda activate

从而复用用户登录 fat2 后已经生效的 conda 初始化配置。


十七、v9.4 inflation factor 卡住修复
-----------------------------------

v9.3 使用：

    bash -ic "bash run_inflation_factor.sh ..."

会让一个交互式 shell 在 SecureCRT/PTTY 环境中长时间占用当前会话，
可能导致网页控制端不断发送状态标记，而 inflation 计算没有真正开始。

v9.4 改为：

1. 仅用一个最长 30 秒的交互式 bash 查询一次：

       conda activate
       command -v python

2. 得到 conda 激活后真正使用的 Python 路径。

3. 检查该 Python 能 import numpy 和 scipy。

4. 返回普通非交互 workflow，用该 Python 逐 trait 运行：

       get_genomicInflationFactor.py

不再让整个 inflation 计算运行在 bash -ic 中。

inflationFactor 目录位置不是项目根目录，而是：

    <project>/EMMAx.Result/hIBS/<archive>/inflationFactor/

例如：

    <project>/EMMAx.Result/hIBS/GPallInd_LocYear_blue/inflationFactor/


十八、v9.5 inflationFactor 输出目录与 conda base
------------------------------------------------

膨胀系数输出目录改为当前项目根目录：

    <project>/inflationFactor/

例如：

    /data9/home/yzhao/GWAS_IRGSP1.0/GPallIndmiss20_UAV_549lines_5years_Re-fitting_Blue/inflationFactor/

每个性状输出：

    <project>/inflationFactor/<trait>.inflationFactor

运行膨胀系数之前，在 fat2 的交互式 shell 中明确执行：

    conda activate base

然后读取该 base 环境中的 Python 路径：

    command -v python

确认 Python 能加载 numpy 和 scipy 后，再返回普通 workflow 批量运行
get_genomicInflationFactor.py。

这样不会让整个 inflation 计算长期占用交互式 bash，但可以保证实际使用
conda base 环境中的 Python。


十九、v9.6 inflation 卡住修复 + inflation-only
----------------------------------------------

1. 不再使用 bash -ic。
   直接利用登录环境中导出的 CONDA_EXE：

       eval "$("$CONDA_EXE" shell.bash hook)"
       conda activate base

   然后用 base 环境 Python 运行 get_genomicInflationFactor.py。

2. 新增 workflow_v2.sh 的 inflation stage。
   如果此前 GWAS 已经生成并归档了所有 .ps 文件，不需要重跑
   prepare / qsub / plot / peak SNP。

   例如 BLUE 结果归档为：

       EMMAx.Result/hIBS/GPallInd_LocYear_blue/

   可直接运行：

       bash workflow_v2.sh inflation \
         GPallInd \
         GPallIndmiss20_UAV_549lines_5years_Re-fitting_Blue \
         GPallInd_LocYear_blue.txt \
         inflation-recovery \
         LocYear-blue

   它读取：

       EMMAx.Result/hIBS/GPallInd_LocYear_blue/traits.txt
       EMMAx.Result/hIBS/GPallInd_LocYear_blue/output/*_chr01_hIBS.ps

   并直接生成：

       <project>/inflationFactor/*.inflationFactor

3. 如果归档 output 或 traits.txt 缺失，会直接报错，不会自动从头重跑。


二十、v9.7 网页 inflation-only 与本地停止
-----------------------------------------

- 历史任务新增“只计算膨胀系数（使用已有 .ps）”按钮。
  它创建新的 inflation-only 恢复任务，直接读取已有归档 output/*.ps，
  不重新运行 prepare / qsub / EMMAX / plot / peak SNP。

- 新增“停止本地监控 / 清除运行状态”按钮。
  只清理 Windows 网页状态，不会自动 qdel PBS，也不会删除服务器结果。

- stop_blue.cmd 现在会自动扫描 runs/*/state.txt，
  将 connecting/running 改为 cancelled，并写入 cancelled.flag。
  不需要再手工 PowerShell 修改 state.txt。

- inflation-only 使用新的 run id，避免旧卡住的 SecureCRT 监控覆盖新任务状态。


二十一、v9.8 inflation Python 固定路径 + VBS 修复
-------------------------------------------------

- 不再使用 conda activate。
- fat2 固定使用：
  /data6/tool/anaconda2-4.1.1/bin/python
- inflation-only 会从服务器 canonical program 目录：
  /data9/home/yzhao/GWAS_IRGSP1.0/000data_prepare/program/get_genomicInflationFactor.py
  复制到当前恢复 run 的 program 目录后运行。
- 正常完整 GWAS/BLUE final 使用 prepare 阶段已经复制到 batch/program 的同一程序。
- 修复 v9.7 inflation-only bridge.vbs 的 VBScript “无效字符”错误。
- 输出仍为：
  <project>/inflationFactor/


二十二、v9.9 单文件原始 GWAS 修复
--------------------------------

问题：
    单个地点/单个年份/单个表型文件时，workflow_blue_v8.sh 会按设计
    委托给原始 workflow_v2.sh。

    v9.8 的 server_blue.py 只在 BLUE 多文件 bundle 中上传 workflow_v2.sh，
    单文件路径提前 return，导致远程报：

        ERROR: missing remote workflow_v2.sh

修复：
    单文件提交时也自动把本地原始 workflow_v2.sh 上传到本次远程 staging：

        /data9/home/yzhao/GWAS_IRGSP1.0/.gwas-web/<runid>/workflow_v2.sh

    然后 workflow_blue_v8.sh 正常委托给原始单表型 GWAS 流程。

多文件 BLUE 逻辑不变：
    仍上传 BLUE-aware patched workflow_v2.sh。

inflation-only、固定 Python、网页停止按钮等 v9.8 功能全部保留。
