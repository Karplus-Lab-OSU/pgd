***************
Search Workflow
***************

PGD Search workflow is designed to take advantage of all of the features provided by Django: Models, Forms, and QuerySets. Search Models and Forms are interchangeable via methods provided by PGD.

------
Models
------

.. .. image:: search_models.png

The Search class closely mimics a protein segment. It is composed of a three classes: Search, SearchResidue, and SearchCode. By creating Search as django models they can be stored and retrieved from the database.

^^^^^^
Search
^^^^^^

Search model contains all fields present in a protein.

^^^^^^^^^^
SearchCode
^^^^^^^^^^

SearchCode is a list of protein IDs to include in the search

^^^^^^^^^^^^^
SearchResidue
^^^^^^^^^^^^^

SearchResidue mimics a Residue object. It contains all fields that a

    * Fields are all strings and can use the PGD Query Syntax
    * Type fields such as AA_Type and SS_Type are stored as an integer with choices encoded in binary to conserve space.

-----
Forms
-----

SearchForms mimic Search Models exactly.

-----------
Conversions
-----------

Search and SearchForm are interchangeable.

    * Use **function name** to convert a **Search** to a **!SearchForm**
    * Use **function name** to convert a **!SearchForm** to a **Search**

Search can generate Django QuerySets

    Use **search.queryset()**

--------
Workflow
--------

.. .. image:: query_lifecycle.png

When submitting the search form, if the search is successful a **Search** object will be created and stored in the user's session. The **Search** can be retrieved and converted back into a **SearchForm** to edit the search.

.. .. image:: query_customization.png

This **Search** object becomes the basis for rendering all other pages. PGD is taking advantage of two features of querysets:

    * **Lazy-execution of SQL Queries** - the query is not performed until the results are requested from the object.
    * **Querysets can be further customized** - the query can be refined further to filter the fields returned, statistical calculations, etc.

These features allow PGD apply view specific logic without recreating the **base query** for every **view**. The **base query** contains the set of records to operate on. The **view** applies its specific logic and calculations.

^^^^^^^^
Examples
^^^^^^^^

    * Plot page only requires 3 properties to perform its calculation
    * Statistics only needs to return Averages, Standard Deviation and other statistics instead of the residue values.

