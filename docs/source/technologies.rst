**********************************************
Technologies Used By Protein Geometry Database
**********************************************

==================
Macro technologies
==================

These technologies are the large building blocks about which the
software is designed.

Updating this software is usually a long, involved process, so no
attempt at documenting this process will be included here.

----------
Django 1.6
----------

Homepage
    https://www.djangoproject.com/
Latest Release
    https://github.com/django/django/releases

Django is a high-level Python Web framework that encourages rapid
development and clean, pragmatic design. Django provides some of the
following features that are useful to PGD:

    * Object Relational Mapping (ORM) - maps python classes to SQL
      tables, including a comprehensive query engine.
    * Form API for validation of user input
    * Templating system for layout of presentation

---------
MySQL 5.5
---------

Homepage
    https://www.mysql.com/
Latest Release
    http://dev.mysql.com/downloads/

MySQL is "the world's most popular open source database".

==================
Micro technologies
==================

These technologies are just as important but are smaller and usually
easier to update and/or replace.  The Javascript tools may be
updated by using the following process (using jquery as an example):

#. Check out a copy of the PGD from Github and create a new branch.

   ::

      $ git clone https://github.com/osuosl/pgd.git
      $ (cd pgd && git checkout develop && git checkout -b update-jquery)

#. Retrieve the latest release.

   ::

      $ wget https://github.com/jquery/jquery/archive/1.11.3.tar.gz


#. Extract the `jquery.min.js` and `jquery.min.map` files and put them
   in the correct location.

   .. note: Not all packages include a map file!

   ::

      $ tar zxf 1.11.3.tar.gz
      $ cp jquery-1.11.3/dist/jquery.min.* pgd/pgd_core/static/js/

   If the filenames have changed, perhaps because the version number
   is included in the filename: 

   #. use `git grep` on the PGD source tree to find all references to
      the old version.  Replace these references with references to
      the new version.  Use `git add` to stage these files for a
      commit.

   #. Use `git rm` to remove the old file, and `git add` to add the
      new file.

#. Spin up an instance, perform a search, and look at the plot to make
   sure the text still looks right.

#. Use `git commit` to commit your changes to your local repository.
   Your commit message should include the new version as well as any
   changes or differences noted during testing.

#. Use `git push` to push your changes to Github.  Finally, submit a
   pull request based on this branch.

------------
jQuery 1.4.1
------------

Homepage
    http://jquery.com
Latest Release
    https://github.com/jquery/jquery/releases

jQuery is a fast and concise JavaScript Library that simplifies HTML
document traversing, event handling, animating, and Ajax interactions
for rapid web development. jQuery is used extensively in the front end
to provide a "web 2.0" experience with dynamically updating pages.

---------------------
jQuery qTip 1.0.0 RC3
---------------------

Homepage
    http://qtip2.com/
Latest Release
    http://qtip2.com/download

qTip is a jQuery plugin which provides extensive support for almost
every kind of tooltip imaginable.

At this time, it is best to retrieve this product from the downloads
page via their "build your own" tool.

* Select the stable build, not the nightly.
* Select all included styles.
* Select all included features, except IE6 support.
* Do not select any additional libraries.

Extract the "jquery.qtip.min.js" file from the downloaded archive and
follow the instructions above.

-------------
Raphael 2.1.4
-------------

Homepage
    http://raphaeljs.com
Latest Release
    https://github.com/DmitryBaranovskiy/raphael/releases

Raphaël is a small JavaScript library that should simplify your work
with vector graphics on the web. Raphael is cross-browser supported
and is used to render graphs within PGD.

==================
Other technologies
==================

These technologies are used to create data files such as fonts which
are used by PGD.

---------------
DejaVu 400 Font
---------------

The DejaVu font included in PGD was generated with `cufon`_.

To rebuild the font:

#. Download the latest version of the font from `Sourceforge`_.
#. Visit the `cufon font generator`_ site.
#. Upload the DejaVuSans files that correspond to the choices.
#. Select "Basic Latin" and "Greek and Coptic".
#. Set the units-per-em to 360.
#. Set the receiving function to "Raphael.registerFont".
#. Acknowledge and accept terms, then hit the "Let's do this!" button.

These instructions are untested!

.. _cufon: https://github.com/sorccu/cufon
.. _Sourceforge: http://sourceforge.net/projects/dejavu/?source=typ_redirect
.. _yui compressor: https://github.com/yui/yuicompressor
.. _cufon font generator: http://cufon.shoqolate.com/generate/


============
Soon to come
============

These technologies are soon to be added to PGD in the hopes of
simplifying things or making them work better.

----------------
jQuery UI 1.11.4
----------------

Homepage
    http://jqueryui.com/
Latest Releases
    http://jqueryui.com/download/

This plugin provides a comprehensive set of user interface widgets,
including autocomplete.  It will replace the existing autocomplete
Javascript and CSS.  It may also replace the tooltips.

Select the `Stable` link under `Quick downloads` on the
latest-releases page.  In addition to copying the `jquery-ui.min.js`
into place, also copy the `jquery-ui.min.css` file into the
appropriate directory for CSS files.


--------------------
jQuery Browser 0.0.7
--------------------

Homepage
    https://github.com/gabceb/jquery-browser-plugin
Latest Release
    https://github.com/gabceb/jquery-browser-plugin/releases

jQuery stopped supporting browser detection in 1.9, so the
functionality is available via a separate plugin.

Alternatively, this functionality may be replaced with feature checks
instead of browser checks.  Have to wait and see!


