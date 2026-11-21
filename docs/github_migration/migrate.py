#!/usr/bin/env python3
"""
This script was used to migrate Torch from Bitbucket to GitHub on June 17, 2026.
It was executed by the machine account torch-sf-assistant to clearly separate
which information was migrated.

Issues and pull requests were downloaded beforehand (now stored in ./files/).

Closed pull requests were uploaded to GitHub as clearly labeled (closed) issues
to preserve documentation of past discussions.

This script was heavily inspired by jeffwidman's migrate.py [1], but was
tailored to better suit our needs.

[1] https://github.com/jeffwidman/bitbucket-issue-migration/blob/master/migrate.py
"""
import os
import sys
import json
import time
import urllib.request
import urllib.error
import random
import re
from collections import defaultdict

GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")
GITHUB_REPO = os.getenv("GITHUB_REPO")
STATE_FILE = os.getenv("STATE_FILE", "migration_state.json")
API_REQUEST_LIMIT = int(os.getenv("API_REQUEST_LIMIT", "60"))
DRY_RUN=True

def load_state():
    if os.path.exists(STATE_FILE):
        try:
            with open(STATE_FILE, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            print(f"Warning: failed to load state: {e}. Starting fresh.")
    return {"issues": {}, "pullrequests": {}}

def save_state(state):
    try:
        with open(STATE_FILE, "w", encoding="utf-8") as f:
            json.dump(state, f, indent=2)
    except Exception as e:
        print(f"Error saving state: {e}", file=sys.stderr)


def make_github_request(url, data=None, token=None, method=None):
    """Helper to send requests to GitHub API using standard library urllib."""
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/vnd.github+json",
        "User-Agent": "Bitbucket-To-GitHub-Migrator",
        "X-GitHub-Api-Version": "2022-11-28"
    }
    
    req_data = None
    if data is not None:
        req_data = json.dumps(data).encode("utf-8")
        headers["Content-Type"] = "application/json"

    # Sleep to respect rate limit
    if API_REQUEST_LIMIT > 0:
        sleep_time = (3600.0 / API_REQUEST_LIMIT) * 1.05
        # Add a tiny random jitter to make it look less bot-like
        jitter = random.uniform(-0.05 * sleep_time, 0.05 * sleep_time)
        actual_sleep = max(0.1, sleep_time + jitter)
        print(f"[Rate Limiter] Sleeping for {actual_sleep:.2f} seconds...")
        time.sleep(actual_sleep)

    req = urllib.request.Request(url, data=req_data, headers=headers, method=method)
    
    try:
        with urllib.request.urlopen(req) as response:
            return response.status, json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        err_msg = e.read().decode("utf-8")
        print(f"HTTP Error {e.code} for {url}: {err_msg}", file=sys.stderr)
        raise

def format_date(date_str):
    if not date_str or not isinstance(date_str, str):
        return "unknown"
    if "T" in date_str:
        return date_str.split("T")[0]
    return date_str

def clean_mentions(text, user_map=None):
    if not text:
        return ""
    
    def replace_id_mention(match):
        account_id = match.group(1)
        if user_map and account_id in user_map:
            return user_map[account_id]
        return account_id
    
    text = re.sub(r'@\{([^}]+)\}', replace_id_mention, text)
    
    text = re.sub(r'(?<!\w)@([a-zA-Z0-9_\-]+)', r'\1', text)

    return text

def format_issue_body(issue, user_map=None):
    reporter = issue.get('reporter')
    display_name = reporter.get('display_name') if isinstance(reporter, dict) else reporter
    created_on = format_date(issue.get('created_on', 'unknown'))
    original_id = issue.get('id')
    status = issue.get('status')
    
    header = (
        f"> [!NOTE]\n"
        f"> **Migrated from Bitbucket**\n"
        f"> - **Original Bitbucket Issue:** `#{original_id}`\n"
        f"> - **Original Author:** {display_name}\n"
        f"> - **Original Date:** {created_on}\n"
        f"> - **Original Status:** {status}\n\n"
    )
    return header + clean_mentions(issue.get('content') or "", user_map)

def format_comment_body(comment, user_map=None):
    user = comment.get('user')
    display_name = user.get('display_name') if isinstance(user, dict) else user
    created_on = format_date(comment.get('created_on', 'unknown'))
    
    header = (
        f"> [!NOTE]\n"
        f"> **Original Comment Author:** {display_name}\n"
        f"> - **Original Date:** {created_on}\n\n"
    )
    return header + clean_mentions(comment.get('content') or "", user_map)

def format_pullrequest_body(pr, user_map=None):

    pr_id = pr.get("id")
    title = pr.get("title")
    state = pr.get("state")
    description = pr.get("description") or ""

    author = pr.get("author")
    display_name = (
        author.get("display_name")
        if isinstance(author, dict)
        else author
    )

    created_on = format_date(pr.get("created_on"))
    updated_on = format_date(pr.get("updated_on"))

    source_branch = pr.get("source_branch", "unknown")
    target_branch = pr.get("target_branch", "unknown")

    header = (
        f"> [!NOTE]\n"
        f"> **Migrated from Bitbucket Pull Request**\n"
        f"> - **Original Bitbucket PR:** `#{pr_id}`\n"
        f"> - **Title:** {title}\n"
        f"> - **Author:** {display_name}\n"
        f"> - **Created:** {created_on}\n"
        f"> - **Updated:** {updated_on}\n"
        f"> - **State:** {state}\n"
        f"> - **Source Branch:** `{source_branch}`\n"
        f"> - **Target Branch:** `{target_branch}`\n\n"
    )

    body = clean_mentions(description, user_map)

    if not body.strip():
        body = "_No description provided._"

    return header + body

def migrate_issues(db_path = "issues/db-2.0.json"):
    
    if not os.path.exists(db_path):
        print(f"Error: {db_path} not found.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {db_path}...")
    with open(db_path, "r", encoding="utf-8") as f:
        db = json.load(f)

    issues = sorted(db.get("issues", []), key=lambda x: x["id"])
    comments = db.get("comments", [])
    
    comments_by_issue = defaultdict(list)
    for c in comments:
        comments_by_issue[c["issue"]].append(c)

    # Sort comments by creation time or ID
    for issue_id in comments_by_issue:
        comments_by_issue[issue_id].sort(key=lambda c: c.get("created_on") or str(c.get("id")))

    print(f"Loaded {len(issues)} issues and {len(comments)} comments.")

    user_map = {}
    def collect_user(u):
        if isinstance(u, dict) and "account_id" in u and "display_name" in u:
            user_map[u["account_id"]] = u["display_name"]

    for issue in issues:
        collect_user(issue.get("reporter"))
        collect_user(issue.get("assignee"))
        for w in issue.get("watchers", []):
            collect_user(w)
        for v in issue.get("voters", []):
            collect_user(v)
    for c in comments:
        collect_user(c.get("user"))

    if DRY_RUN:
        print("\n=== DRY RUN MODE ENABLED ===")

    # Map closed/resolved statuses
    # Open: "new", "open", "on hold"
    # Closed: "resolved", "closed", "invalid", "duplicate", "wontfix"
    closed_statuses = {"resolved", "closed", "invalid", "duplicate", "wontfix"}

    state = load_state()

    for issue in issues:
        issue_id = issue["id"]
        title = issue["title"]
        status = issue["status"]
        body = format_issue_body(issue, user_map)
        
        # Check if already migrated
        if str(issue_id) in state["issues"]:
            print(f"Issue #{issue_id} already migrated (GitHub Issue #{state['issues'][str(issue_id)]}). Skipping.")
            continue

        print(f"\nProcessing Issue #{issue_id}: {title} (Status: {status})")
        
        issue_payload = {
            "title": title,
            "body": body
        }

        if DRY_RUN:
            print(f"[Dry Run] Would POST issue to https://api.github.com/repos/{GITHUB_REPO}/issues")
            print(f"  Payload: {json.dumps(issue_payload, indent=2)}")
            gh_issue_number = 999  # Dummy for dry run
        else:
            url = f"https://api.github.com/repos/{GITHUB_REPO}/issues"
            print(f"Creating issue on GitHub...")
            status_code, resp_data = make_github_request(url, data=issue_payload, token=GITHUB_TOKEN)
            gh_issue_number = resp_data["number"]
            print(f"Created GitHub Issue #{gh_issue_number}")
            
            # Save state immediately
            state["issues"][str(issue_id)] = gh_issue_number
            save_state(state)

            # If the original status was closed/resolved, close the issue on GitHub
            if status in closed_statuses:
                print(f"Original status is '{status}'. Closing GitHub issue #{gh_issue_number}...")
                patch_url = f"https://api.github.com/repos/{GITHUB_REPO}/issues/{gh_issue_number}"
                make_github_request(patch_url, data={"state": "closed"}, token=GITHUB_TOKEN, method="PATCH")
                print(f"Closed GitHub Issue #{gh_issue_number}")

        # Migrate comments for this issue
        issue_comments = comments_by_issue.get(issue_id, [])
        for c in issue_comments:
            comment_body = format_comment_body(c, user_map)
            comment_payload = {"body": comment_body}
            
            if DRY_RUN:
                print(f"[Dry Run] Would POST comment for Issue #{issue_id}")
                print(f"  Payload: {json.dumps(comment_payload, indent=2)}")
            else:
                comment_url = f"https://api.github.com/repos/{GITHUB_REPO}/issues/{gh_issue_number}/comments"
                print(f"Creating comment for Issue #{gh_issue_number}...")
                make_github_request(comment_url, data=comment_payload, token=GITHUB_TOKEN)

def migrate_pullrequests(db_path="pullrequests.json"):

    if not os.path.exists(db_path):
        print(f"Error: {db_path} not found.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {db_path}...")
    with open(db_path, "r", encoding="utf-8") as f:
        prs = json.load(f)

    # Sort for deterministic migration
    prs = sorted(prs, key=lambda x: x["id"])

    print(f"Loaded {len(prs)} pull requests.")

    # Build a map of account_id -> display_name to clean up mentions
    user_map = {}
    for pr in prs:
        for comment in pr.get("comments", []):
            user_id = comment.get("user_id")
            user_name = comment.get("user")
            if user_id and user_name:
                user_map[user_id] = user_name

    if DRY_RUN:
        print("\n=== DRY RUN MODE ENABLED ===")

    state = load_state()
    closed_states = {"MERGED", "DECLINED"}

    for pr in prs:
        pr_id = pr["id"]
        title = pr["title"]
        state_val = pr["state"]

        # Check if already migrated
        if str(pr_id) in state["pullrequests"]:
            print(f"PR #{pr_id} already migrated (GitHub Issue #{state['pullrequests'][str(pr_id)]}). Skipping.")
            continue

        print(f"\nProcessing PR #{pr_id}: {title} ({state_val})")

        body = format_pullrequest_body(pr, user_map)

        issue_payload = {
            "title": f"[Archived PR #{pr_id}] {title}",
            "body": body,
            "labels": ["archived-pr"],
        }

        if DRY_RUN:
            print(f"[Dry Run] Would POST PR-as-issue to GitHub")
            print(json.dumps(issue_payload, indent=2))
            gh_issue_number = 888  # Dummy for dry run
            
            for c in pr.get("comments", []):
                comment_body = format_comment_body(c, user_map)
                print(f"  [Dry Run] Would POST comment for PR #{pr_id}")
                print(f"    Payload: {json.dumps({'body': comment_body}, indent=2)}")

        else:
            url = f"https://api.github.com/repos/{GITHUB_REPO}/issues"

            print("Creating GitHub issue for PR...")
            status_code, resp_data = make_github_request(
                url,
                data=issue_payload,
                token=GITHUB_TOKEN
            )

            gh_issue_number = resp_data["number"]
            print(f"Created GitHub Issue #{gh_issue_number}")

            # Save state immediately
            state["pullrequests"][str(pr_id)] = gh_issue_number
            save_state(state)

            # Migrate comments for this PR
            pr_comments = pr.get("comments", [])
            # Sort comments by creation time or ID
            pr_comments.sort(key=lambda c: c.get("created_on") or str(c.get("id")))
            for c in pr_comments:
                comment_body = format_comment_body(c, user_map)
                comment_payload = {"body": comment_body}
                comment_url = f"https://api.github.com/repos/{GITHUB_REPO}/issues/{gh_issue_number}/comments"
                print(f"Creating comment for Issue #{gh_issue_number}...")
                make_github_request(comment_url, data=comment_payload, token=GITHUB_TOKEN)

            if state_val in closed_states:
                print(f"Closing archived PR issue #{gh_issue_number}...")
                patch_url = f"https://api.github.com/repos/{GITHUB_REPO}/issues/{gh_issue_number}"
                make_github_request(
                    patch_url,
                    data={"state": "closed"},
                    token=GITHUB_TOKEN,
                    method="PATCH"
                ) 

if __name__ == "__main__":

    migrate_issues()    
    migrate_pullrequests()  

    print("\nMigration finished successfully!")
