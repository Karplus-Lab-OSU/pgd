****************
Release Workflow
****************

------------------------------------------------------------------------
Announce the upcoming release to the PGD list four weeks before release.
------------------------------------------------------------------------

Include the following information:

    * version number to be released in x.y.z format
    * rough list of features and bugfixes expected to be included
    * schedule of events (feature freeze, release branch, etc.)

--------------------------------------------------------------
Impose a feature freeze on develop three weeks before release.
--------------------------------------------------------------

Features are no longer permitted to be merged into the develop branch, only bugs.

Release engineer checks each resolved ticket to confirm that the branch was indeed merged into the develop branch.

--------------------------------------------------
Start the release branch two weeks before release.
--------------------------------------------------

Create a branch from develop named 'release/x.y.z' using the version number mentioned in the release announcement.

Log into the dev site and check out the release branch there.

Bugfixes can only be made from and returned to this branch at this time.

--------------------------------------------------
Freeze the release branch one week before release.
--------------------------------------------------

Log into the staging site and check out the release branch there.

Only emergency fixes allowed at this point!

-----------------------------------------------------------------------------
Release the software, close tickets and unfreeze develop on the release date.
-----------------------------------------------------------------------------

Merge the release branch back into master and develop branches.

Log into the production site and check out the master branch there.

All resolved tickets should be closed at this time.

Any existing branches should be rebased from develop before development continues.

Features may now be merged back into develop at this time.
