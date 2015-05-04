#ifndef _HS_SFM_UNIT_TEST_SFM_PIPELINE_REAL_DATA_GENERATOR_HPP_
#define _HS_SFM_UNIT_TEST_SFM_PIPELINE_REAL_DATA_GENERATOR_HPP_

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iterator>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"
#include "hs_sfm/fundamental/linear_8_points_ransac_refiner.hpp"

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _Scalar>
class RealDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;

private:
  const static int feature_dimension_ = 128;

public:
  static Err LoadBundlerOutFile(const std::string& out_path,
                                size_t image_width, size_t image_height,
                                KeysetContainer& keysets,
                                MatchContainer& matches)
  {
    typedef EIGEN_VECTOR(Scalar, 2) Key;
    typedef EIGEN_STD_MAP(size_t, Key) ShuffledKeyset;
    typedef EIGEN_STD_VECTOR(ShuffledKeyset) ShuffledKeysetContainer;

    std::ifstream out_file(out_path.c_str());
    if (!out_file)
    {
      return -1;
    }

    std::string line;
    std::getline(out_file, line);
    std::getline(out_file, line);
    std::stringstream ss(line);
    size_t number_of_images;
    size_t number_of_points;
    ss>>number_of_images>>number_of_points;
    for (size_t i = 0; i < number_of_images * 5; i++)
    {
      std::getline(out_file, line);
    }
    TrackContainer tracks_shuffle(number_of_points);
    ShuffledKeysetContainer keysets_shuffle(number_of_images);
    for (size_t i = 0; i < number_of_points; i++)
    {
      std::getline(out_file, line);
      std::getline(out_file, line);
      std::getline(out_file, line);
      ss.clear();
      ss.str(line);
      size_t number_of_views;
      ss>>number_of_views;
      for (size_t j = 0; j < number_of_views; j++)
      {
         size_t image_id;
         size_t key_id;
         Key key;
         ss>>image_id>>key_id>>key[0]>>key[1];
         key[0] = key[0] + image_width * Scalar(0.5);
         key[1] = image_height * Scalar(0.5) - key[1]/* - Scalar(1)*/;
         if (keysets_shuffle[image_id].find(key_id) ==
             keysets_shuffle[image_id].end())
         {
           keysets_shuffle[image_id][key_id] = key;
         }
         tracks_shuffle[i].push_back(std::make_pair(image_id, key_id));
      }
    }

    keysets.resize(number_of_images);
    for (size_t i = 0; i < number_of_images; i++)
    {
      auto itr_key = keysets_shuffle[i].begin();
      auto itr_key_end = keysets_shuffle[i].end();
      for (; itr_key != itr_key_end; ++itr_key)
      {
        keysets[i].AddKey(itr_key->second);
      }
    }

    for (size_t i = 0; i < number_of_points; i++)
    {
      for (size_t j = 0; j < tracks_shuffle[i].size(); j++)
      {
        size_t image_id = tracks_shuffle[i][j].first;
        size_t key_id = tracks_shuffle[i][j].second;
        auto itr_key_begin = keysets_shuffle[image_id].begin();
        auto itr_key = keysets_shuffle[image_id].find(key_id);
        tracks_shuffle[i][j].second =
          size_t(std::distance(itr_key_begin, itr_key));
      }
    }
    MatchesTracksConvertor convertor;
    MatchContainer matches_initial;
    convertor(tracks_shuffle, matches_initial);

    return FilterMatches(keysets, matches_initial, matches);
  }

  static Err LoadGCPs(const std::string& gcp_path,
                      PointContainer& gcps,
                      TrackContainer& tracks_gcp,
                      KeysetContainer& keysets_gcp,
                      PointContainer& check_points,
                      TrackContainer& tracks_check_point,
                      KeysetContainer& keysets_check_point)
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
    Point mean = Point::Zero();
    for (size_t i = 0; i < tracks_gcp_loose.size(); i++)
    {
      if (tracks_gcp_loose[i].size() > 1)
      {
        mean += gcps_loose[i];
        gcps.push_back(gcps_loose[i]);
        tracks_gcp.push_back(tracks_gcp_loose[i]);
      }
    }
    for (size_t i = 0; i < tracks_check_point_loose.size(); i++)
    {
      if (tracks_check_point_loose[i].size() > 1)
      {
        mean += check_points_loose[i];
        check_points.push_back(check_points_loose[i]);
        tracks_check_point.push_back(tracks_check_point_loose[i]);
      }
    }
    mean /= Scalar(gcps.size() + check_points.size());
    for (size_t i = 0; i < gcps.size(); i++)
    {
      gcps[i] -= mean;
    }

    for (size_t i = 0; i < check_points.size(); i++)
    {
      check_points[i] -= mean;
    }

    return 0;
  }

  static Err FilterMatches(const KeysetContainer& keysets,
                           const MatchContainer& matches_initial,
                           MatchContainer& matches_filtered)
  {
    typedef hs::sfm::fundamental::Linear8PointsRansacRefiner<Scalar> Refiner;
    typedef typename Refiner::Key Key;
    typedef typename Refiner::KeyPair RefinerKeyPair;
    typedef typename Refiner::KeyPairContainer RefinerKeyPairContainer;
    typedef typename Refiner::IndexSet IndexSet;

    auto itr_key_pairs = matches_initial.begin();
    auto itr_key_pairs_end = matches_initial.end();
    double distance_threshold = 32.0;
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
              refiner_key_pairs_refined, inlier_indices, 2048);
      std::cout<<"inlier_indices.size():"<<inlier_indices.size()<<"\n";
      if (inlier_indices.size() > 100)
      {
        hs::sfm::KeyPairContainer key_pairs_refined;
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
};

}
}
}

#endif
