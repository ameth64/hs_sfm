#ifndef _HS_SFM_UNIT_TEST_INCREMENTAL_REAL_DATA_GENERATOR_HPP_
#define _HS_SFM_UNIT_TEST_INCREMENTAL_REAL_DATA_GENERATOR_HPP_

#include <string>
#include <vector>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "hs_image_io/whole_io/image_io.hpp"

#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/fundamental/linear_8_points_ransac_refiner.hpp"

#if 1
#include "hs_sfm/sfm_file_io/matches_saver.hpp"
#endif

#define NO_DETECTION 0
#if NO_DETECTION
#include "hs_sfm/sfm_file_io/keyset_loader.hpp"
#endif


namespace hs
{
namespace sfm
{
namespace incremental
{

template <typename _Scalar, typename _ImageDimension>
class RealDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;
  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;

private:
  const static int feature_dimension_ = 128;

public:
  Err operator() (const std::vector<std::string>& image_paths,
                  const std::string& gcp_path,
                  KeysetContainer& keysets,
                  MatchContainer& matches,
                  PointContainer& gcps,
                  TrackContainer& tracks_gcp,
                  KeysetContainer& keysets_gcp) const
  {
    std::vector<std::string> key_paths;
    std::vector<std::string> descriptor_paths;
    if (DetectFeature(image_paths, key_paths, descriptor_paths, keysets) != 0)
    {
      return -1;
    }

    MatchContainer matches_initial;
    if (MatchFeatures(key_paths, descriptor_paths, matches_initial) != 0)
    {
      return -1;
    }

    if (FilterMatches(keysets, matches_initial, matches) != 0)
    {
      return -1;
    }

    if (LoadGCPs(gcp_path, gcps, tracks_gcp, keysets_gcp) != 0)
    {
      return -1;
    }
#if 1
    std::cout.setf(std::ios::fixed);
    std::cout<<std::setprecision(6);
    for (size_t i = 0; i < gcps.size(); i++)
    {
      std::cout<<"gcps["<<i<<"]:"
               <<gcps[i][0]<<" "<<gcps[i][1]<<" "<<gcps[i][2]<<"\n";
      size_t number_of_views = tracks_gcp[i].size();
      std::cout<<"views:"<<number_of_views<<"\n";
      for (size_t j = 0; j < number_of_views; j++)
      {
        size_t image_id = tracks_gcp[i][j].first;
        size_t key_id = tracks_gcp[i][j].second;
        std::cout<<"image:"<<image_id<<" ";
        std::cout<<"x:"<<keysets_gcp[image_id][key_id][0]<<" "
                 <<"y:"<<keysets_gcp[image_id][key_id][1]<<"\n";
      }
    }

#endif
    return 0;
  }

private:
  Err DetectFeature(const std::vector<std::string>& image_paths,
                    std::vector<std::string>& key_paths,
                    std::vector<std::string>& descriptor_paths,
                    KeysetContainer& keysets) const
  {
    key_paths.clear();
    descriptor_paths.clear();
    keysets.clear();
    if (image_paths.empty()) return -1;
    //找出所有影像文件路径的公共部分。
    auto itr_image_path = image_paths.begin();
    auto itr_image_path_end = image_paths.end();
    std::string base_dir = *itr_image_path;
    UniformSeperators(base_dir);
    ++itr_image_path;
    for (; itr_image_path != itr_image_path_end; ++itr_image_path)
    {
      std::string image_path = *itr_image_path;
      UniformSeperators(image_path);
      std::string::size_type pos = 0;
      while (1)
      {
        std::string::size_type pos_found = image_path.find('/', pos);
        if (pos_found == std::string::npos ||
            pos_found != base_dir.find('/', pos) ||
            image_path.substr(0, pos_found) != base_dir.substr(0, pos_found))
        {
          base_dir = base_dir.substr(0, pos);
          break;
        }
        pos = pos_found + 1;
      }
    }

    //创建需要的目录，获取key文件路径以及descriptor文件路径
    itr_image_path = image_paths.begin();
    for (; itr_image_path != itr_image_path_end; ++itr_image_path)
    {
      std::string image_path = *itr_image_path;
      UniformSeperators(image_path);
      std::string file_path =
        image_path.substr(base_dir.size(),
                          image_path.rfind('.') - base_dir.size());
      std::string directory =
        "intermediate/" +
        file_path.substr(0, file_path.rfind('/') + 1);
      boost::filesystem::create_directories(
        boost::filesystem::path(directory));
      std::string key_path = "intermediate/" + file_path + ".key";
      std::string descriptor_path = "intermediate/" + file_path + ".dsc";
      key_paths.push_back(key_path);
      descriptor_paths.push_back(descriptor_path);
    }

#if !NO_DETECTION
    //检测特征点
    cv::SIFT sift_(6000);
    itr_image_path = image_paths.begin();
    auto itr_key_path = key_paths.begin();
    auto itr_key_path_end = key_paths.end();
    auto itr_descriptor_path = descriptor_paths.begin();
    auto itr_descriptor_path_end = descriptor_paths.end();
    hs::imgio::whole::ImageIO image_io;
    for (; itr_image_path != itr_image_path_end;
         ++itr_image_path, ++itr_key_path, ++itr_descriptor_path)
    {
      //imread doesn't work!And I don't know why!
      //cv::Mat image = cv::imread(*itr_image_path, cv::IMREAD_GRAYSCALE);
      std::cout<<"Detecting:"<<*itr_image_path<<"\n";
      hs::imgio::whole::ImageData image_data;
      if (image_io.LoadImage(*itr_image_path, image_data) != 0) return -1;
      cv::Mat image(image_data.height(), image_data.width(), CV_8UC3,
                    image_data.GetBuffer());
      cv::Mat image_gray;
      cv::cvtColor(image, image_gray, cv::COLOR_RGB2GRAY);
      cv::Mat mask;
      std::vector<cv::KeyPoint> keys;
      cv::Mat descriptors;
      sift_(image_gray, mask, keys, descriptors);

      size_t number_of_keys = keys.size();
      Keyset keyset(number_of_keys);
      std::ofstream key_file(itr_key_path->c_str(), std::ios::out);
      if (!key_file) return -1;
      key_file.setf(std::ios::fixed);
      key_file<<std::setprecision(6);
      key_file<<number_of_keys<<"\n";
      for (size_t i = 0; i < number_of_keys; i++)
      {
        key_file<<keys[i].pt.x<<" "<<keys[i].pt.y<<"\n";
        keyset[i] << keys[i].pt.x,
                     keys[i].pt.y;
      }
      keysets.push_back(keyset);
      key_file.close();

      std::ofstream descriptor_file(itr_descriptor_path->c_str(),
                                    std::ios::out | std::ios::binary);
      if (!descriptor_file) return -1;
      descriptor_file.write((const char*)(descriptors.data),
                            number_of_keys * sizeof(float) * 128);
      descriptor_file.close();
    }
#else
    hs::sfm::fileio::KeysetLoader<Scalar> loader;
    keysets.resize(image_paths.size());
    for (size_t i = 0; i < image_paths.size(); i++)
    {
      loader(key_paths[i], keysets[i]);
    }
#endif

    return 0;
  }

  Err MatchFeatures(const std::vector<std::string>& key_paths,
                    const std::vector<std::string>& descriptor_paths,
                    MatchContainer& matches) const
  {
    matches.clear();
    size_t number_of_images = key_paths.size();
    if (descriptor_paths.size() != number_of_images) return -1;
    for (size_t i = 0; i < number_of_images; i++)
    {
      std::cout<<"Matching train image "<<i<<"\n";
      cv::Mat descriptors_index = LoadDescriptors(key_paths[i],
                                                  descriptor_paths[i]);
      cv::flann::Index index(descriptors_index,
                             cv::flann::KDTreeIndexParams(4));
      for (size_t j = 0; j < i; j++)
      {
        std::cout<<"Matching query image "<<j<<"\n";
        cv::Mat descriptors_match = LoadDescriptors(key_paths[j],
                                                    descriptor_paths[j]);
        cv::Mat indices(descriptors_match.rows, 2, CV_32SC1);
        cv::Mat distances(descriptors_match.rows, 2, CV_32FC1);
        index.knnSearch(descriptors_match, indices, distances,
                        2, cv::flann::SearchParams(128));

        const float match_threshold = 0.6f;

        ImagePair image_pair(i, j);
        KeyPairContainer key_pairs;
        for (size_t k = 0; k < descriptors_match.rows; k++)
        {
          if (distances.at<float>(k, 0) <
              distances.at<float>(k, 1) * match_threshold)
          {
            //Found Match!
            key_pairs.push_back(KeyPair(size_t(indices.at<int>(k, 0)), k));
          }
        }
        if (key_pairs.size() > 16)
        {
          matches[image_pair] = key_pairs;
        }
      }
    }

#if 1
    hs::sfm::fileio::MatchesSaver saver;
    saver("test_matches.txt", matches);
#endif
    return 0;
  }

  Err FilterMatches(const KeysetContainer& keysets,
                    const MatchContainer& matches_initial,
                    MatchContainer& matches_filtered) const
  {
    typedef hs::sfm::fundamental::Linear8PointsRansacRefiner<Scalar> Refiner;
    typedef typename Refiner::Key Key;
    typedef typename Refiner::KeyPair RefinerKeyPair;
    typedef typename Refiner::KeyPairContainer RefinerKeyPairContainer;
    typedef typename Refiner::IndexSet IndexSet;

    auto itr_key_pairs = matches_initial.begin();
    auto itr_key_pairs_end = matches_initial.end();
    Scalar distance_threshold = Scalar(32);
    Refiner refiner;
    for (; itr_key_pairs != itr_key_pairs_end; ++itr_key_pairs)
    {
      RefinerKeyPairContainer refiner_key_pairs_initial;
      RefinerKeyPairContainer refiner_key_pairs_refined;
      IndexSet inlier_indices;
      auto itr_key_pair = itr_key_pairs->second.begin();
      auto itr_key_pair_end = itr_key_pairs->second.end();
      for (; itr_key_pair != itr_key_pair_end; ++itr_key_pair)
      {
        RefinerKeyPair refiner_key_pair;
        refiner_key_pair.first =
          keysets[itr_key_pairs->first.first][itr_key_pair->first];
        refiner_key_pair.second =
          keysets[itr_key_pairs->first.second][itr_key_pair->second];
        refiner_key_pairs_initial.push_back(refiner_key_pair);
      }
      refiner(refiner_key_pairs_initial, distance_threshold,
              refiner_key_pairs_refined, inlier_indices);
      if (inlier_indices.size() > 100)
      {
        KeyPairContainer key_pairs_refined;
        for (size_t i = 0; i < inlier_indices.size(); i++)
        {
          key_pairs_refined.push_back(
            itr_key_pairs->second[inlier_indices[i]]);
        }
        matches_filtered[itr_key_pairs->first] = key_pairs_refined;
      }
    }
    return 0;
  }

  Err LoadGCPs(const std::string& gcp_path,
               PointContainer& gcps,
               TrackContainer& tracks_gcp,
               KeysetContainer& keysets_gcp) const
  {
    using namespace boost::property_tree;
    typedef EIGEN_VECTOR(Scalar, 2) Key;
    typedef EIGEN_STD_VECTOR(Key) KeyContainer;

    std::ifstream xml_file(gcp_path);
    if (!xml_file) return -1;
    ptree pt_root;
    try
    {
      read_xml(xml_file, pt_root);
    }
    catch(xml_parser_error&)
    {
      return -1;
    }

    size_t invalid_uint = std::numeric_limits<size_t>::max();
    Scalar invalid_scalar = std::numeric_limits<Scalar>::max();
    std::string invalid_string = "";

    auto pt_document = pt_root.get_child_optional("document");
    if (!pt_document) return -1;

    auto pt_chunk = pt_document->get_child_optional("chunk");
    if (!pt_chunk) return -1;

    auto pt_cameras = pt_chunk->get_child_optional("cameras");
    if (!pt_cameras) return -1;
    size_t number_of_images = pt_cameras->size();
    keysets_gcp.resize(number_of_images);

    auto pt_markers = pt_chunk->get_child_optional("markers");
    if (!pt_markers) return -1;
    size_t number_of_gcps = pt_markers->size();
    gcps.resize(number_of_gcps);
    tracks_gcp.resize(number_of_gcps);
    auto itr_marker = pt_markers->begin();
    auto itr_marker_end = pt_markers->end();
    for (size_t i = 0; itr_marker != itr_marker_end; ++itr_marker)
    {
      if (itr_marker->first != "marker") continue;
      auto itr_location = itr_marker->second.begin();
      auto itr_location_end = itr_marker->second.end();
      for (; itr_location != itr_location_end; ++itr_location)
      {
        if (itr_location->first != "location") continue;
        size_t camera_id =
          itr_location->second.get(path("<xmlattr>/camera_id", '/'),
                                   invalid_uint);
        if (camera_id == invalid_uint) continue;
        Key key;
        key[0] =
          itr_location->second.get(path("<xmlattr>/x", '/'),
                                   invalid_scalar);
        key[1] =
          itr_location->second.get(path("<xmlattr>/y", '/'),
                                   invalid_scalar);
        keysets_gcp[camera_id].AddKey(key);
        tracks_gcp[i].push_back(
          std::make_pair(camera_id, keysets_gcp[camera_id].size() - 1));
      }
      i++;
    }

    auto pt_ground_control = pt_chunk->get_child_optional("ground_control");
    auto itr_record = pt_ground_control->begin();
    auto itr_record_end = pt_ground_control->end();
    for (; itr_record != itr_record_end; ++itr_record)
    {
      if (itr_record->first != "record") continue;
      size_t marker_id =
        itr_record->second.get(path("<xmlattr>/marker_id", '/'),
                               invalid_uint);
      if (marker_id == invalid_uint) continue;
      gcps[marker_id][0] =
        itr_record->second.get(path("<xmlattr>/x", '/'), invalid_scalar);
      gcps[marker_id][1] =
        itr_record->second.get(path("<xmlattr>/y", '/'), invalid_scalar);
      gcps[marker_id][2] =
        itr_record->second.get(path("<xmlattr>/z", '/'), invalid_scalar);
    }

    return 0;
  }

  static void UniformSeperators(std::string& path)
  {
    auto itr = path.begin();
    auto itr_end = path.end();
    for (; itr != itr_end; ++itr)
    {
      if (*itr == '\\') *itr = '/';
    }
  }

  static cv::Mat LoadDescriptors(const std::string& key_path,
                                 const std::string& descriptor_path)
  {
    cv::Mat descriptors;
    while (1)
    {
      std::ifstream key_file(key_path.c_str(), std::ios::in);
      if(!key_file) break;
      size_t number_of_keys;
      key_file >> number_of_keys;
      key_file.close();
      std::ifstream descriptor_file(descriptor_path,
                                    std::ios::in |std::ios::binary);
      if (!descriptor_file) break;
      descriptors.create(number_of_keys, feature_dimension_, CV_32FC1);
      descriptor_file.read((char*)(descriptors.ptr()),
                           number_of_keys * feature_dimension_ * sizeof(float));
      descriptor_file.close();
      break;
    }
    return descriptors;
  }
};

}
}
}

#endif
