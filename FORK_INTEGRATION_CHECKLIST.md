# Fork and Integration Checklist

## ✅ What You Need to Do

### 1. Fork the Repository on GitHub
- [ ] Go to: https://github.com/KasperThystrup/serovar_detector
- [ ] Click "Fork" button (top-right)
- [ ] Confirm it's created at: https://github.com/Jacques-does-science/serovar_detector

### 2. Push Your Changes to Your Fork
```bash
cd /Users/jacquesolivier/envs/Antigravity/swineotype/third_party/serovar_detector

# Update remote to point to your fork
git remote set-url origin https://github.com/Jacques-does-science/serovar_detector.git

# Commit your changes
git add workflow/scripts/summarize_serovars.R workflow/envs/R.yaml workflow/rules/summarise_serovars.smk
git commit -m "Fix R dependencies and kma_dir bug for swineotype integration

- Add explicit library imports to summarize_serovars.R
- Update R.yaml with all required R packages  
- Fix undefined variable bug: kma_dir -> kma_files
- Add log output to summarise_serovars.smk"

# Push to your fork
git push origin dev
```

### 3. Commit Changes to Main swineotype Repo
```bash
cd /Users/jacquesolivier/envs/Antigravity/swineotype

# Stage the changes
git add swineotype/adapters/app.py scripts/install_swineotype.sh README.md

# Commit
git commit -m "Fix APP serotyping installation and dependencies

- Add peppy dependency to installation script
- Copy serovar_profiles.yaml to config directory for R script
- Update README to use forked serovar_detector with bug fixes
- Install package in editable mode for development"

# Push to your repo
git push origin main
```

## 📋 Summary of Changes

### Changes to `swineotype` repo:
1. **scripts/install_swineotype.sh**: Added `peppy` dependency, changed to editable install
2. **swineotype/adapters/app.py**: Copy serovar_profiles.yaml to config directory
3. **README.md**: Updated to point to your fork of serovar_detector

### Changes to `serovar_detector` fork:
1. **workflow/scripts/summarize_serovars.R**: 
   - Added explicit library imports
   - Fixed `kma_dir` → `kma_files` bug
   - Added debug output
2. **workflow/envs/R.yaml**: Updated with all required R packages
3. **workflow/rules/summarise_serovars.smk**: Added log output

## ✨ After These Steps

Anyone who clones your swineotype repo and follows the README will:
1. Clone your fork of serovar_detector (with fixes)
2. Install all required dependencies including peppy
3. Have a working installation that can run APP serotyping

## 🔄 Optional: Keep Your Fork Updated

To keep your fork in sync with the original repo:
```bash
cd third_party/serovar_detector
git remote add upstream https://github.com/KasperThystrup/serovar_detector.git
git fetch upstream
git merge upstream/dev
```
