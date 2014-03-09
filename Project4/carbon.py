
# In[7]:

x = loadtxt("x.txt")
t = loadtxt("t.txt")
C = loadtxt("C.txt")
C_1 = loadtxt("C_1hour.txt")
C_6 = loadtxt("C_6hour.txt")
C_12 = loadtxt("C_12hour.txt")
C_24 = loadtxt("C_24hour.txt")
C_48 = loadtxt("C_48hour.txt")


# In[10]:

figure()
imshow(C)
colorbar(orientation='horizontal')
xlabel("time $t$")
ylabel("distance $x$")
title("$C(x,t)$")


# Out[10]:

#     <matplotlib.text.Text at 0x109b3ca50>

# image file:

# In[33]:

figure()
plot(x, C_1, 'r', label="$C_{1}(x, 3600)$")
plot(x, C_6, 'g--', label="$C_{6}(x, 21600)$")
plot(x, C_12, 'b-.', label="$C_{12}(x, 43200)$")
plot(x, C_24, color="purple", lw=2, ls=":", label="$C_{24}(x, 86400)$")
plot(x, C_48, color="green", lw=2, ls="-", label="$C_{48}(x, 172800)$")
xlabel("distance $x$ in $mm$")
ylabel("Concentration of Carbon")
title("Carbon Concentrations")
legend(loc=1)


# Out[33]:

#     <matplotlib.legend.Legend at 0x10d9715d0>

# image file:

# In[30]:

C_3mm = loadtxt("C_3mm.txt")
figure()
plot(t, C_3mm, 'b--')
xlabel("time $t$ in seconds")
ylabel("Concentration of Carbon")
title("Carbon Concentrations as a function of time at 3mm distance")


# Out[30]:

#     <matplotlib.text.Text at 0x10e713990>

# image file:

# In[ ]:



