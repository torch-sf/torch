# Repository Migration Summary: Bitbucket to GitHub

This document summarizes the procedure and steps that were taken to migrate the **Torch** repository from Bitbucket to GitHub. The migration was performed by Eric Andersson, Sabrina Appel, Claude Cournoyer-Cloutier, Mordecai-Mark Mac Low and Brooke Polak.

---

## 1. Environment Setup

A dedicated Conda environment was configured to facilitate the repository cleaning and migration process:

```bash
conda create --name torch-migrate
conda activate torch-migrate
conda install pip
pip install git-filter-repo
```

---

## 2. Repository Cleaning and History Rewriting

### Step 2.1: Clone and Backup Creation

A mirror clone of the Bitbucket repository was created to retrieve all branches, tags, and commit history. A backup copy was also preserved:

```bash
git clone --mirror git@bitbucket.org:torch-sf/torch.git
cp -r torch.git torch-backup.git
cd torch.git
```

### Step 2.2: Large and Problematic File Removal

Historical instances of large or unused files (`cube128` and `fake_cube`) were purged across all commits and branches:

```bash
git filter-repo \
  --path cube128 \
  --path ref_data/cube128 \
  --path examples/cube128 \
  --path tests/turbsph/cube128 \
  --path fake_cube \
  --invert-paths
```

*Note: The file cube128 was introduced in the first commit to the repository. This means that all git hashes were rewritten in this step.*

### Step 2.3: Commit Message Standardization

Historical commit messages were updated to standardize references to Bitbucket pull requests:

```bash
git filter-repo --message-callback '
import re
msg = message.decode()
msg = re.sub(r"\(pull request #(\d+)\)", r"(Bitbucket pull request `#\1`)", msg)
return msg.encode()
'
```

### Step 2.4: Branch Renaming

Historical branches were standardized and categorized according to the project's branching strategy using the mapping file located at `./files/branch-rename-map.txt`:

```bash
while read old new; do
  git branch -m "$old" "$new"
done < ./files/branch-rename-map.txt
```

### Step 2.5: Contributor Identity Standardization

Historical author emails and names were consolidated to standardize commit attributions using the mailmap configuration file located at `./files/mailmap.txt`:

```bash
git filter-repo --mailmap ./files/mailmap.txt
```

*(Note: Contributors must add and verify these email addresses in their GitHub accounts to ensure correct commit-to-account linking.)*

---

## 3. Post-Rewrite Verification

The repository structure, history, and author lists were verified prior to uploading:

```bash
# Confirmed large files were removed from the commit history
git log --all -- cube128
git log --all -- fake_cube

# Verified repository size and object count
git count-objects -vH

# Verified that author identities were consolidated correctly
git log --all --format='%an <%ae>' | sort -u
```

---

## 4. GitHub Upload & Issues Migration

### Step 4.1: Pushing to GitHub

An empty repository was created under the `torch-sf` organization on GitHub (without initializing a README, license, or `.gitignore`), and the rewritten repository was pushed as a mirror:

```bash
git remote add github <github-repository-url>
git push --mirror github
```

### Step 4.2: Issues & Pull Requests Migration

Bitbucket issues and closed pull requests (archived as closed GitHub issues) were uploaded to the new repository using the migration script `migrate.py`. The script was executed from the machine account `torch-sf-assistant` so that it would be clear that the migration was performed by the machine account and not a user.

```bash
python3 migrate.py
```

Open pull requests were recreated manually on GitHub.

---

### Step 5.1: GitHub Issues Cleanup

GitHub issues were labeled and organized.

### Step 5.2: Bitbucket Repository Cleanup

The original Bitbucket repository was archived, and its README was updated with a notice directing users to the new GitHub repository.
