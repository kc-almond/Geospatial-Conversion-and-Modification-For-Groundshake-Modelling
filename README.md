# Point to Raster for Groundshake Modelling

A comprehensive GUI-based Python application for converting point shapefiles into continuous raster surfaces, specifically designed for seismic hazard analysis and groundshake modeling. Transform discrete seismic measurement points into interpolated raster datasets for advanced spatial analysis.

## 🌟 Features

### Core Functionality
- **Multi-Parameter Support**: Process PGA, Richter, Wald, and Spectral Acceleration (SA) values simultaneously
- **Flexible Interpolation**: Choose from nearest neighbor, bilinear, or cubic interpolation methods
- **Batch Processing**: Handle multiple shapefiles in a single operation
- **Custom Cell Sizes**: Configure spatial resolution for each parameter independently
- **Clipping Support**: Restrict output to specific geographic boundaries using clipping shapefiles
- **Real-time Logging**: Monitor processing progress with detailed timestamped logs

### User Interface
- **Intuitive GUI**: Built with tkinter for cross-platform compatibility
- **Parameter Configuration**: Individual field mapping and cell size settings for each seismic parameter
- **Progress Tracking**: Visual progress indicators and comprehensive logging

### Data Quality & Performance
- **CRS Handling**: Automatic coordinate reference system management and transformation
- **Compressed Output**: LZW-compressed GeoTIFF format for optimal file sizes
- **NoData Management**: Proper handling of missing values and edge effects

## 🚀 Installation

### Prerequisites
```bash
pip install geopandas pandas rasterio shapely pyproj scipy pyogrio
```

### Required Libraries
- **GeoPandas**: Spatial data manipulation and analysis
- **Rasterio**: Raster I/O and processing
- **SciPy**: Scientific computing and interpolation algorithms
- **Shapely**: Geometric operations
- **PyProj**: Coordinate reference system transformations
- **PyOGRIO**: Fast shapefile reading engine

### Setup
1. Clone the repository:
```bash
git clone https://github.com/yourusername/Point-to-Raster-for-Groundshake-Modelling.git
cd Point-to-Raster-for-Groundshake-Modelling
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Run the application:
```bash
python point_to_raster_app.py
```

## 📊 Supported Data Types

### Input Requirements
- **File Format**: Point shapefiles (.shp)
- **Geometry Type**: Point geometries only
- **Attribute Fields**: Numeric values for seismic parameters

### Default Field data (from PHIVOLCS - ACER Project Omega)
- **PGA**: Peak Ground Acceleration
- **Richter_PE**: Richter Scale Peak Elevation
- **Wald_PEIS**: Wald Peak Estimated Intensity Scale
- **SA(1.0)**: Spectral Acceleration at 1.0 second
- **SA(0.2)**: Spectral Acceleration at 0.2 second

### Output Specifications
- **Format**: GeoTIFF (.tif)
- **Data Type**: Float32 for high precision
- **NoData Value**: NaN for missing areas

## 🔧 Usage Guide

### Basic Workflow
1. **Select Input Files**: Add one or more point shapefiles containing seismic data
2. **Configure Output**: Choose destination folder for raster outputs
3. **Set Parameters**: 
   - Map field names to seismic parameters
   - Configure cell sizes for spatial resolution
   - Select interpolation method
4. **Optional Clipping**: Add boundary shapefile to restrict output extent
5. **Process**: Start conversion and monitor progress in real-time

### Cell Size Configuration
```python
# Common cell size examples:
# High resolution (5m):     0.00004556624
# Medium resolution (30m):  0.00027451108
# Low resolution (9 arc):   0.0025178324
```

### Interpolation Methods
- **Nearest Neighbor**: Preserves original values, creates stepped surfaces
- **Bilinear**: Smooth transitions, ideal for continuous phenomena
- **Cubic**: Very smooth surfaces, may introduce slight artifacts

## 🗺️ Applications

### Seismic Hazard Assessment
- **Earthquake Risk Mapping**: Create continuous risk surfaces from point measurements
- **Shakemap Generation**: Interpolate ground motion parameters for affected areas
- **Liquefaction Analysis**: Model soil liquefaction potential across regions

### Urban Planning & Engineering
- **Building Code Compliance**: Assess seismic design requirements by location
- **Infrastructure Planning**: Evaluate ground motion impacts on critical infrastructure
- **Emergency Response**: Develop response scenarios based on expected ground shaking

### Research & Analysis
- **Comparative Studies**: Analyze different seismic parameters across regions
- **Temporal Analysis**: Track changes in seismic hazard over time
- **Statistical Modeling**: Create input datasets for probabilistic seismic hazard analysis

## 📈 Performance Considerations

### Memory Requirements
- **5m Cell Size**: ~25-30GB RAM (high resolution)
- **30m Cell Size**: ~1-2GB RAM (standard resolution)
- **9 arc-second**: ~100-500MB RAM (regional scale)

### Processing Time
- Depends on:
  - Number of input points
  - Output resolution (cell size)
  - Interpolation method complexity
  - Available system resources

### Optimization Tips
- Use appropriate cell sizes for your analysis scale
- Consider clipping to reduce processing extent
- Process large datasets in batches
- Monitor system resources during processing

## 🔬 Technical Details

### Spatial Interpolation
The application uses SciPy's `griddata` function for bilinear and cubic interpolation, and `cKDTree` for nearest neighbor interpolation. The interpolation process:

1. **Grid Generation**: Creates regular grid based on data bounds and cell size
2. **Value Interpolation**: Estimates values at grid points using selected method
3. **Raster Creation**: Converts interpolated grid to georeferenced raster
4. **Optional Clipping**: Masks output to boundary geometry if specified

### Coordinate System Management
- Automatic CRS detection and preservation
- Coordinate transformation for clipping shapefiles
- Proper handling of geographic vs. projected coordinates

### Error Handling
- Comprehensive input validation
- Graceful handling of missing data
- Detailed error reporting and logging

### Development Setup
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


## 📚 References

PHIVOLCS - Introduction to Earthquake
GeoPandas Documentation
Rasterio Documentation
SciPy Interpolation Guide
USGS Shakemap Documentation

