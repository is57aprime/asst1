# Problem 1
On a Mac M1
* (1.93x speedup from 2 threads)
* (2.34x speedup from 4 threads)
* (3.04x speedup from 6 threads)
* (3.59x speedup from 8 threads)
* (4.06x speedup from 10 threads)
* (5.36x speedup from 20 threads)
* (5.57x speedup from 30 threads)

NumThreads has to be a factor of height (1200)
![Uploading image.png…]()

# Problem 6
* Time measurements conclude that setting assignments >> everything else
* Did not have access to the data file, so used a claude generated script (generate_data_v4.py) to generate a dummy file
* Speedup: 150ms on multithreaded code (5 threads because the dummy file had K=5) and 700ms on regular 
