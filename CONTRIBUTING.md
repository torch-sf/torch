Contributing to Torch
=====================

Issue reports, enhancements, etc are very welcome.  We will do our best to
address issue reports, pull requests, and inquiries, but due to limited
resources we cannot promise replies.

If you have write access to this repository, create a new branch and submit a
pull request from that branch to `main` following the
[Feature Branch Workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow).

If you do not have write access to this repository, you can fork the repo and
submit a pull request following the
[Forking Workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow).


HOWTO submit a pull request
---------------------------

1. Make your changes on a branch separate from `main`.

2. Test that your branch builds and runs successfully on the Netherlands
   SURFsara supercomputer `Snellius`.
   If possible, also build and run on 1-2 more compute clusters, especially
   Graham (Compute Canada), Draco (Drexel University), and Rusty (Flatiron
   Institute).  The respective contacts are Claude, Sean, and Sabrina as of
   Fall 2022.

3. Perform any other tests needed.

4. Open a Bitbucket pull request (PR) to merge your branch into `main`.
   In your description, mention any tests of your code changes.
   And, call attention to changes that may affect other peoples' workflow,
   e.g., if you move important files, change compiler calls, alter Torch logic
   or parameters, etc.

5. Join the next Torch user's group meeting and tell folks about your PR.  You
   may get or solicit feedback, and folks may ask for more time to review.

6. Get >=1 other person to review and approve your PR.

5. After getting approval, merge your PR!


HOWTO coding style
------------------

* Follow the [AMUSE style guide](https://amuse.readthedocs.io/en/latest/reference/style_guide.html)
  if you are modifying the AMUSE interface code.
