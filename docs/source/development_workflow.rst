********************
Development workflow
********************

-------------------------------------------------------------------
Create an issue in the issue tracker if one does not already exist.
-------------------------------------------------------------------

All modifications to the source code should be associated with an issue for tracking purposes.

------------------------------------------------------------------
Create a branch from the develop branch based on the issue number.
------------------------------------------------------------------

The branch should be named based on the issue type and number. For now, all issues are considered bugs for naming purposes. The branch for issue #14109 (a feature) should be named 'bug/14109'. At this time, the issue should be updated with the status 'In Progress'. 

---------------------------------------------------------------------
Write the code, including tests and release notes entry if necessary.
---------------------------------------------------------------------

If the branch is fixing a bug, then a test should be written first if at all possible to confirm the bug exists as well as confirm that the bug fix works.

If the branch makes any user-visible changes, then the release notes (news.html) should be updated to reflect the changes. If you are the first developer to post release notes for the latest development version, add an appropriate header above the existing release notes following the example format -- the release engineer will clean it up if necessary when making the next release. 

---------------------------------
Commit code and update the issue.
---------------------------------

If any billable work is done on an issue, then code should be committed, and the issue should be updated with the number of billable hours spent on the task and a summary of the work that was done. Commits and updates should happen whenever major subtasks are completed and at close of business.

--------------------
Request code review.
--------------------

When the code is complete, the tests if any all pass, and the release notes entry has been added if necessary, mark the issue 'Needs Review' and assign it to another team member for review. That individual will review the code, the tests, the entry, and any other associated changes for accuracy and consistency. If the branch is acceptable, the reviewer will mark the issue 'Resolved'. If the branch is unacceptable, the reviewer will mark the issue 'Needs Work'. In both cases, the reviewer will update the issue with relevant information and reassign the issue back to the original developer. 

-------------------------------------------
Merge changes back into the develop branch.
-------------------------------------------

Once the branch passes review, it should be merged back into develop using 'git merge --no-ff' and then develop should be pushed back into the origin. Once this is complete, the devloper should update the issue to that effect.
