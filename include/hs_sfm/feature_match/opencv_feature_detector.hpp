#ifndef _HS_SFM_FEATURE_MATCH_OPENCV_FEATURE_DETECTOR_HPP_
#define _HS_SFM_FEATURE_MATCH_OPENCV_FEATURE_DETECTOR_HPP_

#include <string>
#include <vector>

#include <opencv2/nonfree/features2d.hpp>

#include "hs_sfm/config/hs_config.hpp"

namespace hs
{
namespace sfm
{
namespace feature_match
{

class HS_EXPORT OpenCVFeatureDetector
{
public:
  OpenCVFeatureDetector(int number_of_features = 0,
                        int number_of_octave_layers = 3,
                        double contrast_threshold = 0.04,
                        double edge_threshold = 10,
                        double sigma = 1.6);

  int Detect(const std::vector<std::string>& image_paths,
             const std::vector<std::string>& key_paths,
             const std::vector<std::string>& descriptor_paths,
             unsigned int pyramid = 0,
             double* complete_ratio = nullptr,
             int* keep_work = nullptr);
private:
  cv::SIFT sift_;
};

}
}
}

#endif
