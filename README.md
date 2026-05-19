# Hydrological Analysis
Generic pipeline that runs standard hydrology algorithms for a given area of interest

## Setup for end user

1. Install Docker.
2. Pull the image from DockerHub.
```
docker pull shubham8625/hydro-pipeline:latest
```
3. Clone the repository.
```
git clone https://github.com/mastershubham/hydro_stack/
cd hydro_stack
```
4. Run the pipeline after mounting the data and the necessary files/folders into the container.
```
docker run -it -w $(pwd) -v $(pwd):$(pwd) <name_of_the_image> bash
python hydrological_analysis.py --shp ./data/masalia_tehsil_boundary.shp --output masalia --grassdb ~/grassdata
```

## Setup for Developers
Some people might want to build their own image locally or maybe make modifications in the image.
First clone/download this Github repository. 

```
git clone https://github.com/mastershubham/hydro_stack
```
Move to the appropriate directory.
```
cd hydro_stack
```

We are using Docker for working in a containerized environment. Ensure that docker is installed on your device and has network connection. 

```
docker build -t <name_of_the_image> .
```
Put up the data somewhere preferably under hydro_stack or some folder under it.

Now run the command as:
```
docker run -it -w $(pwd) -v $(pwd):$(pwd) <name_of_the_image> bash
python hydrological_analysis.py --shp ./data/masalia_tehsil_boundary.shp --output masalia --grassdb ~/grassdata
```

## Visualizing Micro-watershed connectivity

Run docker in the interactive mode.
```
docker run -it -w $(pwd) -v $(pwd):$(pwd) <name_of_the_image> bash
```
Then run the visualization script for micro-watershed connectivity.
```
python visualize_mws_connectivity.py
```
