# ifcopenshell

This script will calculate base quantities for IFC elements based on their geometry. It will generate a list and save it in a CSV file.


## Usage/Examples

You will need IfcOpenShell (I personnally use it with BlenderBIM).

First, you should modify line 6-7 to reference your IFC model and the path to your exported CSV.
```python
model = ifcopenshell.open('C:/path/to/IFC_file.ifc')
exportedCsvFilePathAndName = 'C:/path/to/my_csv.csv'
```

You are now ready to execute script.


## CSV/Excel

Once you CSV file is generated, you can follow instructions in "HOW TO USE" tab in Excel file.



## Screenshots

![App Screenshot](https://github.com/arthur-bellemin/ifcopenshell/blob/main/quantities.png)

