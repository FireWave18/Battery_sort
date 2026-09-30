# Battery_sort

Battery Sort is a command-line C++ utility for sorting and grouping Li-ion / LiFePO4 battery cells before building a pack. You enter the measured capacities of all available cells (and, optionally, their internal resistance), choose the pack configuration (e.g. 3S4P), and the program selects the cells with the smallest spread, groups them so that every parallel group has the closest possible capacity and resistance, calculates pack voltage (min / nominal / max) and estimated capacity (Ah and Wh, based on the weakest group), and optionally saves the results to a timestamped .txt file.

Capacity is the main criterion, resistance is secondary and filters out abnormal cells. The selected cells are not necessarily taken from the middle of the capacity range: the set with the smallest spread is chosen wherever it lies.

Ready-made binaries for Linux and Windows are on the [Releases](https://github.com/FireWave18/Battery_sort/releases/latest) page. To build from source: `g++ -O2 -std=c++17 -o Battery_sort Battery_sort_eng_1_1_0.cpp`
