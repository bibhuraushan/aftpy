# Examples

## Reading a single AFT map file

```python
import aftpy as aft
filepath = "/file/path/AFTmap_20150101_0000.h5"
aftmap = aft.AFTmap(filepath)

# To see the details of the data
info = aftmap.info
print(info)

# To get AFT simulated data
bmap = aftmap.aftmap
print(bmap.shape)

# To plot the map with default color map
aftmap.plot()

# To get the polar field over (60 degree)
pfN, pfS = aftmap.polarfield(latlim=60)

# or

pfN, pfS = aftmap.pfN, aftmap.pfS

# To get complete information in Notebook
print(aftmap)
```

## Reading all AFT map files from a directory and its subdirectories

```python
import aftpy as aft

# Load all files from "/file/root/dir" 
fileroot_dir = "/file/root/dir"
aftdata = aft.AFTmaps(fileroot_dir)

# Generate the list of parametes save it to csv using pnadas DataFrame
outfle = "/Users/bjha/Data/AFT/AFT_parameters_hmifd.csv"
df = aftdata.generate_parameters(outfile=outfle, use_saved=False)

# Get the list of timestamps in ISO format
times = aftdata.timestamps

# Get all the aftmaps as AFTmap object as an ittirator 

aftmaps = aftdata.aftmaps
for aftmap in aftmaps:
    # use aftmap as previous example
    ...

# To calculate Carrington Averaged map for butterfly diagram
crmap = aftdata.cravgmap()


# To convert all maps to png and save in a directory.
aftdata.convert_all()
```



## Other Examples

### Example Usage

---

#### Example 1
```python
import aftpy.aftmap as aft

# Initialize AFTload object
loader = aft.AFTmaps(path="/path/to/aft/maps", filetype="h5")

# Convert all AFT map files to FITS format to '/path/to/converted/maps'
loader.convert_all(convert_to="fits", outpath="/path/to/converted/maps", verbose=True)
```
#### Example 2
```python
import aftpy.getaftdata as aftget
# Initialize AFTdownload object
downloader = aftget.AFTdownload()

# Reload the list of AFT map files
downloader.reload_files()

# Get the list of AFT map files within a specified time range
file_list = downloader.get_list(t0=dt.datetime(2023, 1, 1), t1=dt.datetime(2023, 1, 7))

# Download AFT map files listed in the DataFrame
downloader.download(file_list)
```