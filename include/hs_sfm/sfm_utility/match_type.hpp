#ifndef _HS_SFM_UTILITY_MATCHE_TYPE_HPP_
#define _HS_SFM_UTILITY_MATCHE_TYPE_HPP_

#include <vector>
#include <map>

#include "cereal/types/vector.hpp"
#include "cereal/types/map.hpp"
#include "cereal/types/utility.hpp"

namespace hs
{
namespace sfm
{

/**
 *  表示所有图片对应的与其具有匹配关系的图片id。
 */
typedef std::vector<std::vector<int> > MatchGuide;
/**
 *  表示照片对，first与second分别表示两张照片的id。
 */
typedef std::pair<size_t, size_t> ImagePair;
/**
 *  表示特征点对，first与second分别表示该两个特征点在各自照片中的id。
 */
typedef std::pair<size_t, size_t> KeyPair;
/**
 *  表示多个特征点对。
 */
typedef std::vector<KeyPair> KeyPairContainer;
/**
 *  表示照片对对应的所有匹配点。
 */
typedef std::map<ImagePair, KeyPairContainer> MatchContainer;
/**
 *  表示某三维点在不同照片中的特征点。first表示照片id，second表示特征点id。
 */
typedef std::vector<std::pair<size_t, size_t> > Track;
/**
 *  表示多个Track。
 */
typedef std::vector<Track> TrackContainer;
/**
 *  表示某张照片观测到的所有track的id，first表示track的id，
 *  second表示特征点在track中的id。
 */
typedef std::vector<std::pair<size_t, size_t> > CameraView;
/**
 *  表示多个CameraView。
 */
typedef std::vector<CameraView> CameraViewContainer;

/**
 *  目标索引映射类。
 *
 *  该类表示一组目标（object）与经过映射后的索引（index）之间的关系。
 *  假设一组包含5个元素的目标经过映射后的序列，映射后变成包含4个元素的序列：
 *  1，0，4，此时目标中的第0个元素则对应映射后的序列中的第1个元素，
 *  第1个元素对应映射后的序列中第0个元素，
 *  第5个元素对应映射后的序列中第2个元素，
 *  而第3个以及第4个元素在映射后没有相对应的元素。
 *  此时mapper中的元素为：1,0,invalid,invalid,2。
 */
class ObjectIndexMap
{
public:
  ObjectIndexMap(){}
  ObjectIndexMap(const ObjectIndexMap& other) {mapper_ = other.mapper_;}
  ObjectIndexMap(size_t number_of_objects)
    : mapper_(number_of_objects, invalid_value()) {}
public:

  ObjectIndexMap& operator = (const ObjectIndexMap& other)
  {
    mapper_ = other.mapper_;
    return *this;
  }

  size_t GetMappedId(size_t object_id) const
  {
    return mapper_[object_id];
  }

  void SetObjectId(size_t object_id, size_t mapped_id)
  {
    mapper_[object_id] = mapped_id;
  }

  size_t operator[] (size_t object_id) const
  {
    return mapper_[object_id];
  }

  size_t& operator[] (size_t object_id)
  {
    return mapper_[object_id];
  }

  bool IsValid(size_t object_id) const
  {
    return (mapper_[object_id] != invalid_value());
  }

  void Resize(size_t number_of_objects)
  {
    mapper_.resize(number_of_objects, invalid_value());
  }

  size_t Size() const
  {
    return mapper_.size();
  }

  void AddObject(size_t mapped_id)
  {
    mapper_.push_back(mapped_id);
  }

  bool operator == (const ObjectIndexMap& other) const
  {
    return (mapper_ == other.mapper_);
  }

  template <typename Archive>
  void save(Archive& archive) const
  {
    archive(mapper_);
  }

  template <typename Archive>
  void load(Archive& archive)
  {
    archive(mapper_);
  }

private:
  inline size_t invalid_value() const
  {
    return std::numeric_limits<size_t>::max();
  }

private:
  std::vector<size_t> mapper_;
};

/**
 *  表示一个特征点的可视信息，包括track id、image id和key id。
 *
 */
struct ViewInfo
{
  size_t track_id;
  size_t image_id;
  size_t key_id;
  bool is_blunder;

  bool operator == (const ViewInfo& other) const
  {
    return (track_id == other.track_id &&
      image_id == other.image_id &&
      key_id == other.key_id &&
      is_blunder == other.is_blunder);
  }

  template <typename Archive>
  void serialize(Archive& archive)
  {
    archive(track_id, image_id, key_id, is_blunder);
  }
};

/**
 *  特征点可视信息索引类。提供各种特征点可视信息查询方式。
 */
class ViewInfoIndexer
{
public:
  void Clear()
  {
    views_info_.clear();
    track_image_index_.clear();
    image_key_index_.clear();
  }

  /**
   *  根据所有的track信息创建索引。
   */
  void SetViewInfoByTracks(const hs::sfm::TrackContainer& tracks)
  {
    Clear();
    size_t number_of_tracks = tracks.size();
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = tracks[i].size();
      for (size_t j = 0; j < number_of_views; j++)
      {
        ViewInfo view_info;
        view_info.track_id = i;
        view_info.image_id = tracks[i][j].first;
        view_info.key_id = tracks[i][j].second;
        view_info.is_blunder = false;
        views_info_.push_back(view_info);
        size_t view_id = views_info_.size() - 1;
        track_image_index_[std::make_pair(view_info.track_id,
                                          view_info.image_id)] = view_id;
        image_key_index_[std::make_pair(view_info.image_id,
                                        view_info.key_id)] = view_id;
      }
    }
  }

  /**
   *  根据track id和image id查询特征点可视信息。
   */
  const ViewInfo* GetViewInfoByTrackImage(size_t track_id,
                                          size_t image_id) const
  {
    auto view_itr = track_image_index_.find(std::make_pair(track_id,
                                                            image_id));
    if (view_itr != track_image_index_.end())
    {
      return &views_info_[view_itr->second];
    }
    else
    {
      return nullptr;
    }
  }

  /**
   *  根据track id和image id查询特征点可视信息。
   */
  ViewInfo* GetViewInfoByTrackImage(size_t track_id,
                                    size_t image_id)
  {
    auto view_itr = track_image_index_.find(std::make_pair(track_id,
                                                            image_id));
    if (view_itr != track_image_index_.end())
    {
      return &views_info_[view_itr->second];
    }
    else
    {
      return nullptr;
    }
  }

  /**
   *  根据image id和key id查询特征点可视信息。
   */
  const ViewInfo* GetViewInfoByImageKey(size_t image_id,
                                        size_t key_id) const
  {
    auto view_itr = image_key_index_.find(std::make_pair(image_id, key_id));
    if (view_itr == image_key_index_.end())
    {
      return nullptr;
    }
    else
    {
      return &views_info_[view_itr->second];
    }
  }

  /**
   *  根据image id和key id查询特征点可视信息。
   */
  ViewInfo* GetViewInfoByImageKey(size_t image_id,
                                  size_t key_id)
  {
    auto view_itr = image_key_index_.find(std::make_pair(image_id, key_id));
    if (view_itr == image_key_index_.end())
    {
      return nullptr;
    }
    else
    {
      return &views_info_[view_itr->second];
    }
  }

  bool operator == (const ViewInfoIndexer& other) const
  {
    return (views_info_ == other.views_info_ &&
            track_image_index_ == other.track_image_index_ &&
            image_key_index_ == other.image_key_index_);
  }

  template <typename Archive>
  void serialize(Archive& archive)
  {
    archive(views_info_, track_image_index_, image_key_index_);
  }

private:
  std::vector<ViewInfo> views_info_;
  std::map<std::pair<size_t, size_t>, size_t> track_image_index_;
  std::map<std::pair<size_t, size_t>, size_t> image_key_index_;
};


}
}

#endif
