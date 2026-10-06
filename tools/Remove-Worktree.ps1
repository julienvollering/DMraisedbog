# Remove-Worktree.ps1: safe cleanup of task worktrees under .claude/worktrees ####
#
# Usage (run from the MAIN checkout, not from inside a worktree you are removing):
#   tools\Remove-Worktree.ps1 -List             # status of every worktree, changes nothing
#   tools\Remove-Worktree.ps1 -Name foo, bar    # remove these (if clean and merged)
#   tools\Remove-Worktree.ps1 -AllMerged        # remove every clean, merged worktree
#   tools\Remove-Worktree.ps1 -Leftovers        # clear unregistered folders (links only)
#   add -DryRun to any of the above to see what would happen
#
# SAFETY. A worktree is removed only if (1) nothing is uncommitted, (2) its HEAD is already
# in main, and (3) your shell is not inside it. Each worktree holds `data` and `output`
# junctions to the main checkout's real folders; those are unlinked with a non-recursive
# Directory.Delete BEFORE git removes the folder. Nothing here uses --force, -D, or
# Remove-Item -Recurse.

[CmdletBinding()]
param(
  [string[]]$Name,
  [switch]$AllMerged,
  [switch]$Leftovers,
  [switch]$List,
  [switch]$DryRun
)

$ErrorActionPreference = "Stop"

# Registered worktrees: path + branch, from git's own list. The first entry is the main checkout.
function Get-Worktrees {
  $items = @()
  $cur = $null
  foreach ($line in (git worktree list --porcelain)) {
    if ($line -like "worktree *") {
      if ($cur) { $items += $cur }
      $cur = [pscustomobject]@{ Path = $line.Substring(9); Branch = "" }
    } elseif ($line -like "branch *") {
      $cur.Branch = $line.Substring(7) -replace "^refs/heads/", ""
    }
  }
  if ($cur) { $items += $cur }
  $items
}

function Get-Norm($p) { [System.IO.Path]::GetFullPath($p).TrimEnd("\", "/").ToLower() }

# A real junction or symlink. NOT the ReparsePoint attribute: OneDrive marks ordinary
# synced files and folders as reparse points too (Files On-Demand placeholders), so that
# attribute matches nearly everything in this repository. LinkType is set only for real links.
function Test-IsLink($item) { $item.LinkType -in @("Junction", "SymbolicLink") }

# Links (junctions / symlinks) directly inside a folder.
function Get-Links($dir) {
  Get-ChildItem -LiteralPath $dir -Force |
    Where-Object { Test-IsLink $_ }
}

# Unlink without following: Directory.Delete(path) on a junction removes only the link.
function Remove-Links($dir) {
  foreach ($l in (Get-Links $dir)) {
    if ($DryRun) { Write-Host "    [dry run] unlink $($l.Name)" } else {
      [System.IO.Directory]::Delete($l.FullName, $false)
      Write-Host "    unlinked $($l.Name)"
    }
  }
}

$wts = Get-Worktrees
$main = $wts[0].Path
$root = Get-Norm (Join-Path $main ".claude\worktrees")
$here = Get-Norm (Get-Location).Path

# Task worktrees only: registered, under .claude/worktrees, not the main checkout.
$tasks = @($wts | Select-Object -Skip 1 | Where-Object { (Get-Norm $_.Path).StartsWith($root + "\") })

# Inspect one worktree: returns state and, if it blocks removal, the reason.
function Get-State($wt) {
  if (-not (Test-Path -LiteralPath $wt.Path)) { return [pscustomobject]@{ Ok = $false; Why = "folder missing (run git worktree prune)" } }
  $p = Get-Norm $wt.Path
  if ($here -eq $p -or $here.StartsWith($p + "\")) { return [pscustomobject]@{ Ok = $false; Why = "your shell is inside it" } }
  $dirty = @(git -C $wt.Path status --porcelain).Count
  if ($dirty -gt 0) { return [pscustomobject]@{ Ok = $false; Why = "$dirty uncommitted change(s)" } }
  git -C $wt.Path merge-base --is-ancestor HEAD main 2>$null
  if ($LASTEXITCODE -ne 0) { return [pscustomobject]@{ Ok = $false; Why = "has commits not in main" } }
  [pscustomobject]@{ Ok = $true; Why = "clean and merged" }
}

function Remove-OneWorktree($wt) {
  $leaf = Split-Path $wt.Path -Leaf
  $s = Get-State $wt
  if (-not $s.Ok) { Write-Host "SKIP ${leaf}: $($s.Why)"; return }
  Write-Host "REMOVE $leaf (branch $($wt.Branch))"
  Remove-Links $wt.Path
  if ($DryRun) { Write-Host "    [dry run] git worktree remove; git branch -d"; return }
  git worktree remove $wt.Path
  if ($LASTEXITCODE -ne 0) { Write-Host "    git worktree remove failed; leaving branch alone"; return }
  if ($wt.Branch) { git branch -d $wt.Branch }
  # A locked or half-removed folder can survive; remove it only if truly empty.
  if ((Test-Path -LiteralPath $wt.Path) -and -not (Get-ChildItem -LiteralPath $wt.Path -Force)) {
    Remove-Item -LiteralPath $wt.Path
  }
  if (Test-Path -LiteralPath $wt.Path) { Write-Host "    folder still exists (locked by a process or OneDrive?)" }
}

if ($List -or (-not $Name -and -not $AllMerged -and -not $Leftovers)) {
  foreach ($wt in $tasks) {
    $s = Get-State $wt
    "{0,-28} {1,-6} {2}" -f (Split-Path $wt.Path -Leaf), $(if ($s.Ok) { "OK" } else { "KEEP" }), $s.Why
  }
  # Unregistered folders still on disk.
  $registered = $tasks | ForEach-Object { Get-Norm $_.Path }
  foreach ($d in (Get-ChildItem -LiteralPath (Join-Path $main ".claude\worktrees") -Directory -Force)) {
    if ($registered -notcontains (Get-Norm $d.FullName)) {
      "{0,-28} {1,-6} {2}" -f $d.Name, "LEFT", "not a registered worktree; use -Leftovers"
    }
  }
  return
}

if ($Name) {
  foreach ($n in $Name) {
    $wt = $tasks | Where-Object { (Split-Path $_.Path -Leaf) -eq $n }
    if (-not $wt) { Write-Host "SKIP ${n}: no such registered worktree"; continue }
    Remove-OneWorktree $wt
  }
}

if ($AllMerged) {
  foreach ($wt in $tasks) { Remove-OneWorktree $wt }
}

if ($Leftovers) {
  $registered = $tasks | ForEach-Object { Get-Norm $_.Path }
  foreach ($d in (Get-ChildItem -LiteralPath (Join-Path $main ".claude\worktrees") -Directory -Force)) {
    if ($registered -contains (Get-Norm $d.FullName)) { continue }
    $p = Get-Norm $d.FullName
    if ($here -eq $p -or $here.StartsWith($p + "\")) { Write-Host "SKIP $($d.Name): your shell is inside it"; continue }
    $nonlinks = @(Get-ChildItem -LiteralPath $d.FullName -Force | Where-Object { -not (Test-IsLink $_) })
    if ($nonlinks.Count -gt 0) { Write-Host "SKIP $($d.Name): holds real files ($($nonlinks.Name -join ', ')); inspect by hand"; continue }
    Write-Host "CLEAR $($d.Name) (links only)"
    Remove-Links $d.FullName
    if (-not $DryRun) { Remove-Item -LiteralPath $d.FullName }
  }
}

if (-not $DryRun) { git worktree prune }
