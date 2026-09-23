GWAS v32 - 数据版本前端同步 / runs 稳定性 / 提交错误提示
========================================================

本次只需要替换：
server.py

workflow_v2.sh 继续使用 v31，不需要修改。

根据本次诊断结果
----------------
后端 /api/config 已经正确返回：

GPall    -> GPallmiss20-864lines
GPallInd -> GPallIndmiss20-864lines
GPallJap -> GPallJapmiss20-864lines

因此“数据版本”网页仍为空，不是后端数据源配置错误，而是网页前端保存/
重绘状态没有正确采用后端默认值。

同时磁盘 runs/ 中 meta.json、state.txt 均存在，而且直接调用 snapshot()
能够读取旧 run，因此历史任务“偶尔全部消失”属于本地状态文件并发读取/
API 静默跳过问题，不是运行记录被删除。

v32 修改
---------
1. “数据版本”输入框
   - 优先通过 sampleList 输入框的位置反查数据版本输入框；
   - 也保留“数据版本”标签和“基因型数据源”区域的定位方式；
   - 空值、miss20-v1、旧的自动群体版本都会自动更新为当前群体版本；
   - 自定义版本名称仍保留，不强行改写。

2. 真正提交时再次兜底
   即使网页内部保存状态仍然是空版本，发往：
   /api/preview
   /api/run
   的 JSON 在发送前也会自动补成当前群体的正确 version。

   因此“输入框显示问题”不会再影响真正提交的数据源版本。

3. /api/runs 不再静默丢历史任务
   原来 snapshot() 遇到 OSError/ValueError 会直接 continue，
   整个 run 会从网页历史列表消失。

   v32 改为：
   - bridge 状态文件读取遇到 Windows 短暂占用时，每50ms重试，最多20次；
   - 即使 snapshot 最终失败，也返回该 run 的 meta；
   - 状态显示 attention；
   - detail 显示 Local snapshot error，而不是把 run 隐藏。

4. “开始分析”失败时明确提示
   /api/run 返回错误时，网页右上角会直接显示：
   提交失败：<server.py 返回的具体错误>

   提交成功时显示：
   已创建新任务：<run ID>

   这样不会再出现“点了没反应，下面仍然停留在旧日志”的假象。

保留功能
--------
- v31 population-specific 数据版本和 phenotype 提前处理
- 同名 phenotype 网页确认
- 确认后精确覆盖、不批量清理
- v28 下载日志文件名带群体
- v27 结果目录和 sigSNP 命名
- v26 bridge 状态文件写入重试
- qsub 间隔 1 秒
- Windows 完成/失败通知
- SecureCRT 审计日志

替换
----
1. 双击 stop.cmd
2. 覆盖 server.py
3. 双击 start.cmd
4. 浏览器 Ctrl+F5

检查版本：
Invoke-RestMethod 'http://127.0.0.1:8765/api/config' | Select-Object version

应显示：
32

已完成 server.py Python 语法检查。
