#ifndef _HS_SFM_UNIT_TEST_INCREMENTAL_REAL_SYNTHETIC_DATA_GENERATOR_HPP_
#define _HS_SFM_UNIT_TEST_INCREMENTAL_REAL_SYNTHETIC_DATA_GENERATOR_HPP_

#include <string>
#include <fstream>
#include <sstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/projective_functions.hpp"
#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
{

template <typename _Scalar>
class RealSyntheticDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef size_t ImageDimension;
  typedef int Err;
  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;
  typedef typename Keyset::Key Key;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;
  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef hs::sfm::ImageParams<ImageDimension> Image;
  typedef EIGEN_STD_VECTOR(Image) ImageContainer;

private:
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

public:
  RealSyntheticDataGenerator(
    Scalar outlier_ratio,
    Scalar key_stddev)
    : outlier_ratio_(outlier_ratio)
    , key_stddev_(key_stddev) {}

  Err Generate(const std::string& bundler_out_path,
               const std::string& gcp_path,
               const IntrinsicParamsContainer& intrinsic_params_set_true,
               const std::vector<size_t>& image_intrinsic_map,
               ImageDimension image_width, ImageDimension image_height,
               hs::sfm::MatchContainer& matches,
               PointContainer& points_absolute_true,
               ExtrinsicParamsContainer& extrinsic_params_set_absolute_true,
               KeysetContainer& keysets_noised,
               PointContainer& gcps,
               hs::sfm::TrackContainer& tracks_gcp,
               KeysetContainer& keysets_gcp_noised,
               PointContainer& check_points,
               hs::sfm::TrackContainer& tracks_check_point,
               KeysetContainer& keysets_check_point_noised) const
  {
    Err result = 0;
    while (1)
    {
      KeysetContainer keysets_true;
      result = LoadBundlerOutFile(bundler_out_path,
                                  intrinsic_params_set_true,
                                  image_intrinsic_map,
                                  extrinsic_params_set_absolute_true,
                                  points_absolute_true,
                                  keysets_true,
                                  matches);
      if (result != 0) break;

      Image image;
      image.m_width = image_width;
      image.m_height = image_height;
      ImageContainer images(extrinsic_params_set_absolute_true.size(),
                            image);

      result = GenerateNoisedKeysets(keysets_true,
                                     images,
                                     keysets_noised);
      if (result != 0) break;

      result = LoadGCPs(gcp_path,
                        gcps,
                        tracks_gcp,
                        keysets_gcp_noised,
                        check_points,
                        tracks_check_point,
                        keysets_check_point_noised);
      if (result != 0) break;

      Point mean = Point::Zero();
      for (size_t gcp_id = 0; gcp_id < gcps.size(); gcp_id++)
      {
        mean += gcps[gcp_id];
      }
      mean /= Scalar(gcps.size());

      for (size_t gcp_id = 0; gcp_id < gcps.size(); gcp_id++)
      {
        gcps[gcp_id] -= mean;
      }

      for (size_t check_point_id = 0; check_point_id < check_points.size();
           check_point_id++)
      {
        check_points[check_point_id] -= mean;
      }

      for (size_t point_id = 0; point_id < points_absolute_true.size();
           point_id++)
      {
        points_absolute_true[point_id] -= mean;
      }

      for (size_t extrinsic_id = 0;
           extrinsic_id < extrinsic_params_set_absolute_true.size();
           extrinsic_id++)
      {
        extrinsic_params_set_absolute_true[extrinsic_id].position() -= mean;
      }

      break;
    }

    return result;
  }

  Scalar key_stddev() const
  {
    return key_stddev_;
  }

private:
  Err LoadBundlerOutFile(
    const std::string& bundler_out_path,
    const IntrinsicParamsContainer& intrinsic_params_set_true,
    const std::vector<size_t>& image_intrinsic_map,
    ExtrinsicParamsContainer& extrinsic_params_set_absolute_true,
    PointContainer& points_absolute_true,
    KeysetContainer& keysets_true,
    hs::sfm::MatchContainer& matches) const
  {
    std::ifstream bundler_out_file(bundler_out_path);
    if (!bundler_out_file)
    {
      return -1;
    }

    std::stringstream ss;
    std::string line;
    std::getline(bundler_out_file, line);
    std::getline(bundler_out_file, line);
    ss.str(line);
    size_t number_of_images;
    size_t number_of_points;
    ss>>number_of_images>>number_of_points;

#ifdef GETLINE
#undef GETLINE
#endif
#define GETLINE \
        std::getline(bundler_out_file, line);\
        ss.clear();\
        ss.str(line);

    hs::sfm::ObjectIndexMap image_extrinsic_map(number_of_images);
    for (size_t image_id = 0; image_id < number_of_images; image_id++)
    {
      GETLINE;
      Scalar focal;
      ss>>focal;
      if (focal == 0) continue;
      Matrix33 rotation_matrix;
      GETLINE;
      ss>>rotation_matrix(0, 0)>>rotation_matrix(0, 1)>>rotation_matrix(0, 2);
      GETLINE;
      ss>>rotation_matrix(1, 0)>>rotation_matrix(1, 1)>>rotation_matrix(1, 2);
      GETLINE;
      ss>>rotation_matrix(2, 0)>>rotation_matrix(2, 1)>>rotation_matrix(2, 2);
      GETLINE;
      Point translate;
      ss>>translate[0]>>translate[1]>>translate[2];
      ExtrinsicParams extrinsic_params;
      extrinsic_params.position() = -rotation_matrix.transpose() * translate;
      rotation_matrix.row(1) *= Scalar(-1);
      rotation_matrix.row(2) *= Scalar(-1);
      extrinsic_params.rotation() = rotation_matrix;
      extrinsic_params_set_absolute_true.push_back(extrinsic_params);
      image_extrinsic_map[image_id] =
        extrinsic_params_set_absolute_true.size() - 1;
    }

    hs::sfm::TrackContainer tracks;
    keysets_true.resize(extrinsic_params_set_absolute_true.size());
    for (size_t point_id = 0; point_id < number_of_points; point_id++)
    {
      GETLINE;
      Point point;
      ss>>point[0]>>point[1]>>point[2];
      GETLINE;
      GETLINE;
      size_t number_of_views;
      ss>>number_of_views;
      hs::sfm::Track track;
      for (size_t view_id = 0; view_id < number_of_views; view_id++)
      {
        size_t image_id, key_id;
        Scalar key_x, key_y;
        ss>>image_id>>key_id>>key_x>>key_y;
        key_x = 3000 + key_x;
        key_y = 2000 - key_y;
        if (image_extrinsic_map.IsValid(image_id))
        {
          size_t extrinsic_id = image_extrinsic_map[image_id];
          const ExtrinsicParams& extrinsic_params =
            extrinsic_params_set_absolute_true[extrinsic_id];
          size_t intrinsic_id = image_intrinsic_map[image_id];
          const IntrinsicParams& intrinsic_params =
            intrinsic_params_set_true[intrinsic_id];

          Key key = ProjectiveFunctions<Scalar>::WorldPointProjectToImageKey(
              intrinsic_params, extrinsic_params, point);

          keysets_true[extrinsic_id].AddKey(key);
          track.push_back(
            std::make_pair(image_id, keysets_true[extrinsic_id].size() - 1));
        }
      }

      if (track.size() > 1)
      {
        tracks.push_back(track);
        points_absolute_true.push_back(point);
      }
    }

#undef GETLINE

    hs::sfm::MatchesTracksConvertor convertor;
    convertor(tracks, matches);

    return 0;
  }

  Err GenerateNoisedKeysets(const KeysetContainer& keysets_true,
                            const ImageContainer& images,
                            KeysetContainer& keysets_noised,
                            bool generate_outliers = true) const
  {
    keysets_noised = keysets_true;
    EIGEN_MATRIX(Scalar, 2, 2) covariance_key;
    covariance_key.setIdentity();
    covariance_key *= key_stddev_ * key_stddev_;
    size_t number_of_keysets = keysets_true.size();
    for (size_t i = 0; i < number_of_keysets; i++)
    {
      size_t number_of_keys = keysets_true[i].size();
      for (size_t j = 0; j < number_of_keys; j++)
      {
        Scalar random;
        hs::math::random::UniformRandomVar<Scalar, 1>::Generate(
          0, 1, random);
        if (random < outlier_ratio_ && generate_outliers)
        {
          Key min_key;
          min_key.setZero();
          Key max_key;
          max_key << Scalar(images[i].m_width),
                     Scalar(images[i].m_height);
          hs::math::random::UniformRandomVar<Scalar, 2>::Generate(
            min_key, max_key, keysets_noised[i][j]);
        }
        else
        {
          Key mean = keysets_true[i][j];
          hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
            mean, covariance_key, keysets_noised[i][j]);
        }
      }
    }

    return 0;
  }

  Err LoadGCPs(const std::string& gcp_path,
               PointContainer& gcps,
               hs::sfm::TrackContainer& tracks_gcp,
               KeysetContainer& keysets_gcp,
               PointContainer& check_points,
               hs::sfm::TrackContainer& tracks_check_point,
               KeysetContainer& keysets_check_point) const
  {
    using namespace boost::property_tree;

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
    size_t number_of_images = 0;
    auto itr_cameras_child = pt_cameras->begin();
    auto itr_cameras_child_end = pt_cameras->end();
    for (; itr_cameras_child != itr_cameras_child_end; ++itr_cameras_child)
    {
      if (itr_cameras_child->first == "group")
      {
        auto itr_camera = itr_cameras_child->second.begin();
        auto itr_camera_end = itr_cameras_child->second.end();
        for (; itr_camera != itr_camera_end; ++itr_camera)
        {
          if (itr_camera->first == "camera")
          {
            number_of_images++;
          }
        }
      }
      else if (itr_cameras_child->first == "camera")
      {
        number_of_images++;
      }
    }

    PointContainer gcps_loose;
    TrackContainer tracks_gcp_loose;
    PointContainer check_points_loose;
    TrackContainer tracks_check_point_loose;
    keysets_gcp.resize(number_of_images);
    keysets_check_point.resize(number_of_images);

    auto pt_markers = pt_chunk->get_child_optional("markers");
    if (!pt_markers) return -1;
    size_t number_of_gcps = pt_markers->size();
    auto itr_marker = pt_markers->begin();
    auto itr_marker_end = pt_markers->end();
    std::map<size_t, size_t> marker_gcp_map;
    std::map<size_t, size_t> marker_check_point_map;
    for (; itr_marker != itr_marker_end; ++itr_marker)
    {
      if (itr_marker->first != "marker") continue;
      size_t marker_id = itr_marker->second.get(path("<xmlattr>/id", '/'),
                                                invalid_uint);
      if (marker_id == invalid_uint) continue;
      auto pt_reference = itr_marker->second.get_child_optional("reference");
      if (pt_reference)
      {
        Point point;
        point[0] = pt_reference->get(path("<xmlattr>/x", '/'),
                                     invalid_scalar);
        point[1] = pt_reference->get(path("<xmlattr>/y", '/'),
                                     invalid_scalar);
        point[2] = pt_reference->get(path("<xmlattr>/z", '/'),
                                     invalid_scalar);
        std::string enable = pt_reference->get(path("<xmlattr>/enabled", '/'),
                                               invalid_string);
        if (enable == "true")
        {
          marker_gcp_map[marker_id] = gcps_loose.size();
          gcps_loose.push_back(point);
        }
        else if (enable == "false")
        {
          marker_check_point_map[marker_id] = check_points_loose.size();
          check_points_loose.push_back(point);
        }
      }
    }

    tracks_gcp_loose.resize(gcps_loose.size());
    tracks_check_point_loose.resize(check_points_loose.size());
    auto pt_frames = pt_chunk->get_child_optional("frames");
    if (!pt_frames) return -1;
    auto pt_frame = pt_frames->get_child_optional("frame");
    if (!pt_frame) return -1;
    auto pt_frame_markers = pt_frame->get_child_optional("markers");
    if (!pt_frame_markers) return -1;
    auto itr_frame_marker = pt_frame_markers->begin();
    auto itr_frame_marker_end = pt_frame_markers->end();
    for (; itr_frame_marker != itr_frame_marker_end; ++itr_frame_marker)
    {
      if (itr_frame_marker->first != "marker") continue;
      size_t marker_id =
        itr_frame_marker->second.get(path("<xmlattr>/marker_id", '/'), invalid_uint);
      if (marker_id == invalid_uint) continue;
      auto itr_location = itr_frame_marker->second.begin();
      auto itr_location_end = itr_frame_marker->second.end();
      for (; itr_location != itr_location_end; ++itr_location)
      {
        if (itr_location->first != "location") continue;
        size_t camera_id =
          itr_location->second.get(path("<xmlattr>/camera_id", '/'),
                                   invalid_uint);
        if (camera_id == invalid_uint) continue;
        Key key;
        key[0] = itr_location->second.get(path("<xmlattr>/x", '/'), invalid_scalar);
        key[1] = itr_location->second.get(path("<xmlattr>/y", '/'), invalid_scalar);
        if (marker_gcp_map.find(marker_id) != marker_gcp_map.end())
        {
          size_t gcp_id = marker_gcp_map[marker_id];
          keysets_gcp[camera_id].AddKey(key);
          size_t key_id = keysets_gcp[camera_id].size() - 1;
          tracks_gcp_loose[gcp_id].push_back(
            std::make_pair(camera_id, key_id));
        }
        else if(marker_check_point_map.find(marker_id) !=
                marker_check_point_map.end())
        {
          size_t check_point_id = marker_check_point_map[marker_id];
          keysets_check_point[camera_id].AddKey(key);
          size_t key_id = keysets_check_point[camera_id].size() - 1;
          tracks_check_point_loose[check_point_id].push_back(
            std::make_pair(camera_id, key_id));
        }
      }
    }

    gcps.clear();
    tracks_gcp.clear();
    check_points.clear();
    tracks_check_point.clear();
    for (size_t i = 0; i < tracks_gcp_loose.size(); i++)
    {
      if (tracks_gcp_loose[i].size() > 1)
      {
        gcps.push_back(gcps_loose[i]);
        tracks_gcp.push_back(tracks_gcp_loose[i]);
      }
    }
    for (size_t i = 0; i < tracks_check_point_loose.size(); i++)
    {
      if (tracks_check_point_loose[i].size() > 1)
      {
        check_points.push_back(check_points_loose[i]);
        tracks_check_point.push_back(tracks_check_point_loose[i]);
      }
    }
    return 0;
  }

private:
  Scalar outlier_ratio_;
  Scalar key_stddev_;
};

}
}
}

#endif
