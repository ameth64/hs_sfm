#include <fstream>
#include <iomanip>

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "hs_sfm/feature_match/opencv_feature_detector.hpp"

namespace hs
{
namespace sfm
{
namespace feature_match
{

OpenCVFeatureDetector::OpenCVFeatureDetector(int number_of_features,
                                             int number_of_octave_layers,
                                             double contrast_threshold,
                                             double edge_threshold,
                                             double sigma)
  : sift_(number_of_features, number_of_features, contrast_threshold,
          edge_threshold, sigma)
{
}

int OpenCVFeatureDetector::Detect(
  const std::vector<std::string>& image_paths,
  const std::vector<std::string>& key_paths,
  const std::vector<std::string>& descriptor_paths,
  unsigned int pyramid,
  double* complete_ratio,
  int* keep_work)
{
  size_t number_of_images = image_paths.size();
  if (number_of_images != key_paths.size()) return -1;
  if (number_of_images != descriptor_paths.size()) return -1;

  double scale = 1.0;
  while (pyramid > 0)
  {
    scale /= 2.0;
    pyramid--;
  }

  for (size_t i = 0; i < number_of_images; i++)
  {
    if (keep_work)
    {
      if (!*keep_work)
      {
        break;
      }
    }

    cv::Mat image = cv::imread(image_paths[i].c_str(), cv::IMREAD_GRAYSCALE);
    cv::Mat resized_image;
    if (pyramid > 0)
    {
      cv::resize(image, resized_image, cv::Size(), scale, scale);
    }
    else
    {
      resized_image = image;
    }
    cv::Mat mask;
    std::vector<cv::KeyPoint> keys;
    cv::Mat descriptors;
    sift_(resized_image, mask, keys, descriptors);

    size_t number_of_keys = keys.size();
    std::ofstream key_file(key_paths[i].c_str(), std::ios::out);
    if (!key_file) return -1;
    key_file.setf(std::ios::fixed);
    key_file<<std::setprecision(6);

    key_file.close();

    std::ofstream descripor_file(descriptor_paths[i].c_str(),
                                 std::ios::out | std::ios::binary);
    if (!descripor_file) return -1;
    descripor_file.write((const char*)(descriptors.data),
                         number_of_keys * sizeof(float) * 128);
    descripor_file.close();
    if (complete_ratio)
    {
      *complete_ratio = double(i + 1) / double(number_of_images);
    }
  }

  return 0;
}

}
}
}
