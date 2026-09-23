$ErrorActionPreference = 'Stop'
$projectDirectory = $PSScriptRoot
$pythonCommand = (Get-Command python.exe).Source
$pythonWindowless = Join-Path (Split-Path $pythonCommand) 'pythonw.exe'
if (-not (Test-Path -LiteralPath $pythonWindowless)) { throw "Missing pythonw.exe: $pythonWindowless" }
$startupDirectory = [Environment]::GetFolderPath('Startup')
$shortcutPath = Join-Path $startupDirectory 'GWAS Local Console.lnk'
$shell = New-Object -ComObject WScript.Shell
$shortcut = $shell.CreateShortcut($shortcutPath)
$shortcut.TargetPath = $pythonWindowless
$shortcut.Arguments = '"' + (Join-Path $projectDirectory 'launcher.py') + '"'
$shortcut.WorkingDirectory = $projectDirectory
$shortcut.WindowStyle = 7
$shortcut.Description = 'GWAS local webpage at http://127.0.0.1:8765/ (no automatic analysis)'
$shortcut.Save()
Start-Process -FilePath $pythonWindowless -ArgumentList $shortcut.Arguments -WorkingDirectory $projectDirectory -WindowStyle Hidden
Write-Output "Installed current-user startup shortcut: $shortcutPath"
Write-Output 'Bookmark: http://127.0.0.1:8765/'
