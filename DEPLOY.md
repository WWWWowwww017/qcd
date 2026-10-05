# Deployment Guide

## Repository Layout

- `index.html`: bilingual research overview, method figures and Application results.
- `code-lab.html`: read-only archived replay workbench.
- `assets/`: logos, four method SVGs and the Application Figure 3 raster.
- `docs/`: `Foundation.pdf` and `Application.pdf`.
- `code/`: public source, four test fixtures, preview JSON, embedded workbench data and `claims.json`.
- `vendor/`: local Plotly, Prism and the small Lucide icon modules used by the page.
- `.github/workflows/`: validation and GitHub Pages workflows.

The site is static. There is no build step, package manager, framework or runtime server requirement.

## GitHub Pages

1. Replace the site URL placeholders before publishing.
2. Put this directory at the repository root and push it to `main`.
3. In **Settings → Pages**, choose **GitHub Actions** as the build source.
4. Wait for `Validate static research site` and then `Deploy static research site` to pass.

### First Release Command

Run this from the site root in PowerShell 7, or use the exact BOM-free implementation below in Windows PowerShell 5.1. Do not replace it with the shorter `Get-Content -Raw` form: Windows PowerShell 5.1 may decode these UTF-8 files as the local ANSI code page and `Set-Content -Encoding utf8` may add a BOM. The replacement URL must end with `/`.

```powershell
$u = 'https://<你的地址>/'   # 结尾必须带 /
$utf8NoBom = New-Object System.Text.UTF8Encoding($false)
Get-ChildItem -Recurse -File -Path . -Include *.html,*.xml,*.txt,*.md |
  Where-Object {
    $_.FullName -notmatch '\\(vendor|code|\.github|\.git)\\' -and
    $_.Name -notin @('DEPLOY.md', 'RELEASE_CHECKLIST.md', 'SCIENTIFIC_REVIEW.md', 'PROJECT_AUDIT.md', 'PRE_DEPLOYMENT_TEST_REPORT.md')
  } |
  ForEach-Object {
    $p = $_.FullName
    $t = [System.IO.File]::ReadAllText($p, [System.Text.Encoding]::UTF8)
    [System.IO.File]::WriteAllText($p, $t.Replace('__SITE_URL__', $u), $utf8NoBom)
  }
rg -n "__SITE_URL__" . --glob '!vendor/**' --glob '!code/**' --glob '!.github/**' --glob '!.git/**' --glob '!DEPLOY.md' --glob '!RELEASE_CHECKLIST.md' --glob '!SCIENTIFIC_REVIEW.md' --glob '!PROJECT_AUDIT.md' --glob '!PRE_DEPLOYMENT_TEST_REPORT.md'
```

The final `rg` command uses the same exclusions and must return no matches in the publishable HTML/XML/TXT files after replacement; the five excluded review documents intentionally retain the literal as documentation. The command is therefore idempotent: a second run makes no further changes. Before the address is known, the repository intentionally retains `__SITE_URL__` in canonical, social, sitemap, robots and 404 metadata. The same command preserves UTF-8 text and writes no BOM; verify the title and language text after running it. The PowerShell 5.1-safe `ReadAllText`/`WriteAllText` form above is required; do not use the shorter `Get-Content -Raw` form.

## Release Attachment

Attach the original `qcd-code-archive.zip` to the versioned GitHub Release. It is not part of the Pages tree and was not re-packed in this round. It contains the larger and historical MATLAB/Python archive, fixtures and saved result files. The public tree contains the lightweight source, preview results and the four fixtures needed by the documented Python tests. The archive's `__MACOSX` and `.DS_Store` entries remain pending author confirmation; do not silently replace the original attachment.

## Scientific Boundary

The homepage figure and four reported values are the π / a₁(1260) results of `Application.pdf`, using the paper's Λ range `[5, 10] GeV²`. The ηc and J/ψ workbench displays are archived replays from the same research programme using Λ = 30 GeV²; they are not reconstructions reported in that paper. The browser does not recompute the numerical results.

## License and Citation

The page is supplied for academic exchange. Cite `Foundation.pdf` and `Application.pdf` for the research claims, and preserve the source/result provenance when redistributing code or replays. Third-party vendor files retain their own notices and licenses.

## Reproducible Release Package

From the site root, the following PowerShell command creates a package outside the site directory. It excludes the audit/review documents, `.github` and `.git`, the Release ZIP, source maps, Python bytecode and the old unreferenced preview JSON while retaining the runtime `code/lab-data.js`, `.nojekyll` and all referenced public files.

```powershell
$root = (Resolve-Path .).Path
$stage = Join-Path $env:TEMP 'qcd-inverse-site-release-stage'
$out = Join-Path (Resolve-Path ..).Path 'qcd-inverse-site-release.zip'
Remove-Item -LiteralPath $stage -Recurse -Force -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Force -Path $stage | Out-Null
$files = Get-ChildItem -LiteralPath $root -Recurse -File -Force | Where-Object {
  $_.FullName -notmatch '[\\/](\.git|\.github|__pycache__)([\\/]|$)' -and
  $_.Extension -notin @('.map', '.pyc', '.pyo') -and
  $_.Name -notin @('PROJECT_AUDIT.md', 'RELEASE_CHECKLIST.md', 'SCIENTIFIC_REVIEW.md', 'PRE_DEPLOYMENT_TEST_REPORT.md', 'qcd-code-archive.zip', 'eta-c-matlab-preview.json')
}
foreach ($file in $files) {
  $relative = $file.FullName.Substring($root.Length + 1)
  $destination = Join-Path $stage $relative
  New-Item -ItemType Directory -Force -Path (Split-Path $destination) | Out-Null
  Copy-Item -LiteralPath $file.FullName -Destination $destination
}
Remove-Item -LiteralPath $out -Force -ErrorAction SilentlyContinue
Add-Type -AssemblyName System.IO.Compression.FileSystem
[System.IO.Compression.ZipFile]::CreateFromDirectory($stage, $out, [System.IO.Compression.CompressionLevel]::Optimal, $false)
Get-ChildItem -LiteralPath $stage -Recurse -File | Measure-Object -Property Length -Sum
Get-Item -LiteralPath $out | Select-Object FullName, Length
Remove-Item -LiteralPath $stage -Recurse -Force
```

The release package is a static upload artifact; the full research archive is deliberately kept as a separate Release attachment. `PROJECT_AUDIT.md`, `RELEASE_CHECKLIST.md`, `SCIENTIFIC_REVIEW.md` and `PRE_DEPLOYMENT_TEST_REPORT.md` are internal review documents and are not staged for Pages or included in the release package.
