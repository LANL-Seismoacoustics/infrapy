#!/usr/bin/env python
# coding: utf-8

# # Establishing a Database Connection

# A series of steps illustrating how to build a database, establish a database connection, test the database connection and pull data utilizing Pisces.

# ## Building a Database

# Data in miniseed or sac formats can be loaded into a sqlite database using commands from pisces. In this tutorial, we will build a sqlite database from a series of SAC files contained in test/data/cli

# In[1]:


get_ipython().system(' pisces sac2db -h')


# In[2]:


get_ipython().system('pisces sac2db sqlite:///test.sqlite ../data/cli/*.SAC')


# As infrapy is an array processing tool, after your sqlite database is created, you will need to update the REFSTA for each array using update_refsta.py

# In[5]:


get_ipython().run_line_magic('run', '../../scripts/update_refsta.py sqlite:///test.sqlite FSU')
get_ipython().run_line_magic('run', '../../scripts/update_refsta.py sqlite:///test.sqlite HWU')
get_ipython().run_line_magic('run', '../../scripts/update_refsta.py sqlite:///test.sqlite LCM')
get_ipython().run_line_magic('run', '../../scripts/update_refsta.py sqlite:///test.sqlite PSU')
get_ipython().run_line_magic('run', '../../scripts/update_refsta.py sqlite:///test.sqlite WMU')


# ## Querying the Database

# We can now connect to our database and practice querying

# In[6]:


import pisces as ps

session = ps.db_connect(backend='sqlite', instance='test.sqlite')


# In[7]:


session


# The tables used in this tutorial are all standard CSS3.0 tables, we define these tables in order to query the database we have created.

# In[8]:


from pisces.tables.css3 import Site, Wfdisc, Origin


# In[9]:


q = session.query(Site).all()


# In[10]:


q


# In[11]:


refsta = session.query(Site).filter(Site.refsta == 'FSU').all()


# In[12]:


refsta


# ### Retrieve and plot waveforms from a query

# In[13]:


from pisces import wfdisc2trace, read_waveform


# In[14]:


get_ipython().run_line_magic('matplotlib', 'notebook')
for wf in session.query(Wfdisc).filter(Wfdisc.chan == 'EDF').limit(4):
    tr = wfdisc2trace(wf)
    data = read_waveform(wf.dir + '/' + wf.dfile, wf.datatype, wf.foff, wf.nsamp)
    tr.plot()


# In[ ]:





# In[ ]:




