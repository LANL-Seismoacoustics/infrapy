.. _schema:

================
Schema
================

The purpose of this document is to define the schema used for the
operation of the infrasound analysis tool, infrapy. The tables described
by this document extend the CSS3.0 or KB core schema to include
information required for the operation of infrapy. This document is
divided into three sections, the first being this introduction. Section
two defines eight new, infrasonic data processing-specific database
tables. Both internal (ORACLE) and external formats for the attributes
are defined, along with a short description of each attribute. Section
three of the document shows the relationships between the different
tables by using entity-relationship diagrams.

This schema is a work in progress and may be updated as development of
infrapy continues.

_______________________
Table Descriptions
_______________________

This section describes the logical structure of each table used in the
Infrapy software package. The name of the table is first, followed by a
description of the purpose and use of the table. Below the description
is a listing of the columns, in the order which they are defined in the
tables. The storage column gives the actual ORACLE datatype for the
column in question. The external format and character positions columns
are provided for the convenience of database users who wish to transfer
data between the ORACLE database tables and flat files.

---------------
Conventions
---------------

The following conventions are used, following Carr et al (2002):

+-----------------------+-----------------------+-----------------------+
| **Element**           | **Appearance**        | **Example**           |
+=======================+=======================+=======================+
| Database table        | Bold                  | **arrival**           |
+-----------------------+-----------------------+-----------------------+
| Database columns      | Italic                | *sta*                 |
+-----------------------+-----------------------+-----------------------+
| Database table and    | Bold.italic           | **arrival**\ *.sta*   |
| column when written   |                       |                       |
| in the dot notation   |                       |                       |
+-----------------------+-----------------------+-----------------------+
| Value of a key or     | Courier font          | arid                  |
| component of a key    |                       |                       |
+-----------------------+-----------------------+-----------------------+

--------------------------------------------------
Table Definitions: Infrapy-Specific Tables
--------------------------------------------------
Table descriptions can be found in the API section of the documentation.
