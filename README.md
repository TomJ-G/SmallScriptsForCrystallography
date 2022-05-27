# Small Scripts For Crystallography

This ia a collection of simple and short scripts I used in everyday crystallographic work. 
At the moment I'll upload all scripts succesively in a separate files and later, I'll compile them in one library.
Feel free to use the code that is posted here.

## Description of scripts
Coming soon!

### fcfoplot.py
This is helpful for comparing multiple Fc-Fo values. By default the logarithmic scale is shown, to emphasize the low value region at which biggest disagreements between Fc and Fo are usually seen.

**How to use it?**

In this example I want to compare two plots. I assign paths to csv files to variables and then pass them to function. Also I set limits to 0.01 - 200.

    P = r"C:\structure\data1.csv"
    P2 = r"C:\structure\data2.csv"
    FcFplot([P,P2],0.01,200)
    
![image](https://user-images.githubusercontent.com/59794882/170627703-6f6fdd96-0d3d-476b-9448-371e6de605ec.png)

