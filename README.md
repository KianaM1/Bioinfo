# Solution to Assignment 1

Creating a New Directory

```
mkdir class
cd class
touch readme.md
```

# Samtools Version
```
apt install samtools
samtools --version
```

Samtools version is 1.19.2

# Making a Nested Directory Structure (in the class directory)
```
mkdir -p coding/files
```

# To Create Files in Different Directories
```
touch coding/files/file1.txt
touch coding/file2.txt
```

# Absolute Path
```
root/class/coding/files/file1.txt
root/class/coding/file2.txt
```

# Relative Path
```
sudo cat ../files/file1.txt
cd ..
sudo cat ../coding/file2.txt
```
