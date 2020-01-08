# Raster

## Download 

Download the folder :

```bash
git clone https://github.com/gwendalp/Raster.git
```


## Compiling

Let's get started ! :page_facing_up:
You are now able to compile the project
You can choose your *.txt file under this format :

```
# latitude    longitude  elevation
48.29762887 -004.41737340 14.334
48.29762971 -004.41735997 14.379
48.29763698 -004.41738809 14.452
...
```

By default, you will find in build/ Brest.txt. Then you chose the number of pixels over the X-Axis, in the example, it's 150 pixels. 

```bash
chmod +x buid.sh
./build.sh
cd build/
./create_raster Brest.txt 150
```

Now you are able to, open the generated binary *.ppm.

## Documentation

You can enjoy a documentation by typing :

```bash
cd build/doc/html
firefox index.html &
```

## Tests

You can test functions :

```bash
cd build
make test
```
Please do not hesitate to add your unit tests.

## Authors

* **Antonin Lizé** :tennis:
* **Gwendal Prise** :ocean:
