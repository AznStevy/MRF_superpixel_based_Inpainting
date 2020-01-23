# MRF_superpixel_based_Inpainting

A Markov Random Field superpixel based algorithm for image inpainting. Method based on Cheng et al. Paper can be found [here](https://www.sciencedirect.com/science/article/pii/S0165168418302883).

Below is a demonstration of how this code works. The original image (left) is masked with a single color (middle) and the algorithm outputs the inpainted version. Metrics about this method can be found in the writeup [here](https://github.com/AznStevy/MRF_superpixel_based_Inpainting/blob/master/Paper/Xu_Stephen_MRF_Inpainting_Final_Paper.pdf).

![Original](https://github.com/AznStevy/MRF_superpixel_based_Inpainting/blob/master/images/fibers.png)
![Original masked](https://github.com/AznStevy/MRF_superpixel_based_Inpainting/blob/master/output/fibers_masked.png)
![Output using Superpixels](https://github.com/AznStevy/MRF_superpixel_based_Inpainting/blob/master/output/fibers_blended_sp.png)
![Original](https://github.com/AznStevy/MRF_superpixel_based_Inpainting/blob/master/images/pebbles.png)
![Original masked](https://github.com/AznStevy/MRF_superpixel_based_Inpainting/blob/master/output/pebbles_masked.png)
![Output using Superpixels](https://github.com/AznStevy/MRF_superpixel_based_Inpainting/blob/master/output/pebbles_blended_sp.png)

# How to Run
* Install Matlab 2017b or higher
* run `mrf_image_inpainting_w_direction_structure_dist_analysis_color.m`
