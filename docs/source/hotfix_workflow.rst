***************
Hotfix Workflow
***************

------------------------------------------------------------
One or more show-stopping bugs are detected between releases
------------------------------------------------------------

An example of a show-stopping bug would be something like #14613 where an unexpected side-effect of a bugfix caused incomplete results to be returned on unrelated searches.

--------------------------------------------
Create a hotfix branch off the master branch
--------------------------------------------

The hotfix branch name will take the form hotfix/x.y.z+1

If the current release is 1.2 the hotfix branch is then hotfix/1.2.1 
If the current release is 3.7.2 the hotfix branch is then hotfix/3.7.3

-----------------------------------------
Create bug branches off the hotfix branch
-----------------------------------------

For each bug that must be fixed, an issue is created and a bug branch named after that issue like all other bugs is created from the hotfix branch instead of the develop branch like all other bugs.

-------------------------------------------------------------
Merge resolved bugs back to hotfix branch and test on staging
-------------------------------------------------------------

As each individual bug is resolved, its branch is merged back into the hotfix branch. The hotfix branch can then be updated on the staging server for testing.

-----------------------------------------------------------------------------
When all bugs are resolved and merged, merge hotfix branch into master branch
-----------------------------------------------------------------------------

Do not forget to increment the version, create a new tag, and update the news page with all bug fixes!

--------------------------------------------
Pull master on production and restart Apache
--------------------------------------------

This should make the new version accessible to the user community.

-------------------------------------
Merge hotfix branch back into develop
-------------------------------------

Once production is back up and running, take the time to merge the hotfix branch back into develop.
