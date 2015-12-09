#ifndef _HS_SFM_SFM_UTILITY_MATCHES_TRACKS_CONVERTOR_HPP_
#define _HS_SFM_SFM_UTILITY_MATCHES_TRACKS_CONVERTOR_HPP_

#include <map>
#include <set>
#include <algorithm>
#include <queue>

#include "hs_sfm/sfm_utility/match_type.hpp"

namespace hs
{
namespace sfm
{

struct MatchesTracksConvertor
{
  typedef int Err;

  struct KeyPairCompair	//关键点比较functor, 以图像序号为比对
  {
    bool operator () (const KeyPair& key_pair0, const KeyPair& key_pair1) const
    {
      return key_pair0.first < key_pair1.first;
    }
  };

  Err operator() (const MatchContainer& matches,
                  TrackContainer& tracks) const	//该重载将match集合转换为track集合, tracks为输出变量.
  {
#if 0
    typedef std::pair<size_t, size_t> View;
    typedef std::map<View, size_t> ViewTrackIndexer;
    typedef std::vector<std::pair<View, size_t> > ViewTrackContainer;
    tracks.clear();
    MatchContainer::const_iterator itr_image_pair = matches.begin();
    MatchContainer::const_iterator itr_image_pair_end = matches.end();
    ViewTrackIndexer view_track_indexer;
    size_t number_of_tracks = 0;
    ViewTrackContainer error_view_tracks;
    for (; itr_image_pair != itr_image_pair_end; ++itr_image_pair)
    {
      size_t image0_id = itr_image_pair->first.first;
      size_t image1_id = itr_image_pair->first.second;
      KeyPairContainer::const_iterator itr_key_pair =
        itr_image_pair->second.begin();
      KeyPairContainer::const_iterator itr_key_pair_end =
        itr_image_pair->second.end();
      for (; itr_key_pair != itr_key_pair_end; ++itr_key_pair)
      {
        size_t key0_id = itr_key_pair->first;
        size_t key1_id = itr_key_pair->second;
        View view0(image0_id, key0_id);
        View view1(image1_id, key1_id);
        ViewTrackIndexer::iterator itr_view0 = view_track_indexer.find(view0);
        ViewTrackIndexer::iterator itr_view1 = view_track_indexer.find(view1);
        if (itr_view0 != view_track_indexer.end())
        {
          if (itr_view1 != view_track_indexer.end())
          {
            if (itr_view0->second != itr_view1->second)
            {
              size_t track0_id = itr_view0->second;
              size_t track1_id = itr_view1->second;
              error_view_tracks.push_back(std::make_pair(view0, track0_id));
              error_view_tracks.push_back(std::make_pair(view1, track1_id));
            }
          }
          else
          {
            size_t track_id = itr_view0->second;
            view_track_indexer[view1] = track_id;
            tracks[track_id].push_back(view1);
          }
        }
        else
        {
          if (itr_view1 != view_track_indexer.end())
          {
            size_t track_id = itr_view1->second;
            view_track_indexer[view0] = track_id;
            tracks[track_id].push_back(view0);
          }
          else
          {
            view_track_indexer[view0] = number_of_tracks;
            view_track_indexer[view1] = number_of_tracks;
            number_of_tracks++;
            Track track;
            track.push_back(view0);
            track.push_back(view1);
            tracks.push_back(track);
          }
        }
      }
    }

    auto itr_error_view_track = error_view_tracks.begin();
    auto itr_error_view_track_end = error_view_tracks.end();
    std::cout<<"Error tracks!\n";
    for (; itr_error_view_track != itr_error_view_track_end;
         ++itr_error_view_track)
    {
      size_t track_id = itr_error_view_track->second;
      std::cout<<"error track "<<track_id<<"\n";
      for (size_t i = 0; i < tracks[track_id].size(); i++)
      {
        std::cout<<tracks[track_id][i].first<<" "
                 <<tracks[track_id][i].second<<"\n";
      }
      tracks[track_id].clear();
    }

    return 0;
#else
    typedef std::map<size_t, std::set<size_t> > ImageNeighborIndexer;
    typedef std::pair<size_t, size_t> View;

    //Normalize matches.Make it symmetric.Sort every key pairs in image pairs.
    MatchContainer matches_normalized = matches;
    auto itr_image_pair = matches.begin();
    auto itr_image_pair_end = matches.end();
    for (; itr_image_pair != itr_image_pair_end; ++itr_image_pair)	//开始遍历匹配对
    {
      ImagePair image_pair = itr_image_pair->first;
      ImagePair image_pair_reverse = ImagePair(image_pair.second,
                                               image_pair.first);	//对调两幅图像序号
      if (matches_normalized.find(image_pair_reverse) ==
          matches_normalized.end())	//如果对调后的图像序号在原match container中不存在, 则认为该匹配对还不是"对称的", 需补足对称的值
      {
        KeyPairContainer key_pairs_reverse;
        for (size_t i = 0; i < itr_image_pair->second.size(); i++)
        {
          const KeyPairContainer& key_pairs = itr_image_pair->second;
          KeyPair key_pair_reverse = KeyPair(key_pairs[i].second,
                                             key_pairs[i].first);
          key_pairs_reverse.push_back(key_pair_reverse);
        }
        matches_normalized[image_pair_reverse] = key_pairs_reverse;
      }
    }

    auto itr_image_pair_normalized = matches_normalized.begin();
    auto itr_image_pair_normalized_end = matches_normalized.end();
    for (; itr_image_pair_normalized != itr_image_pair_normalized_end;
         ++itr_image_pair_normalized)
    {
      std::sort(itr_image_pair_normalized->second.begin(),
                itr_image_pair_normalized->second.end());
    }//排序

    //Build image pair indexer.
    ImageNeighborIndexer image_neighbor_indexer;	//一个Indexer(索引器)的作用是, 对于任一图像编号, 索引到与之配对的所有图像编号, 存入一个set
    itr_image_pair_normalized = matches_normalized.begin();
    itr_image_pair_normalized_end = matches_normalized.end();
    for (; itr_image_pair_normalized != itr_image_pair_normalized_end;
         ++itr_image_pair_normalized)
    {
      size_t image0_id = itr_image_pair_normalized->first.first;
      size_t image1_id = itr_image_pair_normalized->first.second;	//取出每个匹配对的两个图像ID
      auto itr_image0_neighbor = image_neighbor_indexer.find(image0_id);
      if (itr_image0_neighbor == image_neighbor_indexer.end())	//若第一个ID在索引器中不存在, 则添加之
      {
        std::set<size_t> neighbors;
        neighbors.insert(image1_id);
        image_neighbor_indexer[image0_id] = neighbors;
      }
      else
      {
        itr_image0_neighbor->second.insert(image1_id);
      }
      auto itr_image1_neighbor = image_neighbor_indexer.find(image1_id);	//反向操作
      if (itr_image1_neighbor == image_neighbor_indexer.end())
      {
        std::set<size_t> neighbors;
        neighbors.insert(image0_id);
        image_neighbor_indexer[image1_id] = neighbors;
      }
      else
      {
        itr_image1_neighbor->second.insert(image0_id);
      }
    }

    itr_image_pair_normalized = matches_normalized.begin();
    itr_image_pair_normalized_end = matches_normalized.end();
    std::set<View> marked_views;	//基于两视图的计算, View是个pair<size_t, size_t>,存储图像ID及
    std::set<View> error_views;
    hs::sfm::TrackContainer tracks_loose;
    for (; itr_image_pair_normalized != itr_image_pair_normalized_end;
         ++itr_image_pair_normalized)	//遍历match(匹配对)集合
    {
      auto itr_key_pair = itr_image_pair_normalized->second.begin();	//取出该match中的key_pair向量.
      auto itr_key_pair_end = itr_image_pair_normalized->second.end();
      size_t image0_id = itr_image_pair_normalized->first.first;		//取出该匹配对的键中的两个图像ID
      size_t image1_id = itr_image_pair_normalized->first.second;
      for (; itr_key_pair != itr_key_pair_end; ++itr_key_pair)	//遍历该匹配对的所有点集
      {
        size_t key0_id = itr_key_pair->first;
        size_t key1_id = itr_key_pair->second;
        View view0(image0_id, key0_id);	//构造View对象, 空间点在某一ID的图像中的点视图, 即图像ID及其中某一同名点
        if (marked_views.find(view0) != marked_views.end())	//如果已该view已被marked, 则跳过
        {
          continue;
        }
        marked_views.insert(view0);
        std::set<View> track_set;
        track_set.insert(view0);
        std::queue<View> view_queue;
        view_queue.push(view0);
        std::set<size_t> image_touched;
        image_touched.insert(image0_id);
        while (!view_queue.empty())
        {
          View view = view_queue.front();
          view_queue.pop();
          auto itr_neighbors = image_neighbor_indexer.find(view.first);	//索引到当前图像ID的所有邻近匹配图像
          if (itr_neighbors != image_neighbor_indexer.end())	//若当前图像ID存在匹配的邻近图像
          {
            auto itr_neighbor = itr_neighbors->second.begin();
            auto itr_neighbor_end = itr_neighbors->second.end();
            for (; itr_neighbor != itr_neighbor_end; ++itr_neighbor)	//遍历其所有匹配的邻近图像
            {
              size_t image_id = *itr_neighbor;

              ImagePair image_pair(view.first, image_id);
              auto itr_image_pair =
                matches_normalized.find(ImagePair(view.first, image_id));	//从两个邻近图像ID构造一匹配对, 在match集合中寻找, 后面将要用到该匹配对的key pair集合
              if (itr_image_pair == matches_normalized.end())
              {
                return -1;
              }
              KeyPair key_pair_to_find(view.second, 0);	//以当前图像ID上同名点编号为起点, 构造key_pair. 注意下方KeyPairCompair比对方法仅以pair的第一成员为键, 第二成员可忽略
              auto search_result =
                std::equal_range(itr_image_pair->second.begin(),
                                 itr_image_pair->second.end(),
                                 key_pair_to_find, KeyPairCompair());	//在邻近图像组成的match的key pair向量中查找同值范围.
              if (search_result.first == search_result.second)	//若不存在与key_pair_to_find等值的区间, 则跳过;
              {
                continue;
              }
              if (std::distance(search_result.first, search_result.second) > 1)	//计算两个迭代器之间的距离, 即序号差, 若>1表示存在图像1的同名点在图像2有多个匹配, 无效.
              {
                continue;
              }
              View view_new(image_id, search_result.first->second);	//用邻近图像ID和在match的key_pair向量中寻找到的最近pair的邻对值, 组成邻近图像的View.
              if (track_set.find(view_new) != track_set.end())	//如果该View在track_set中存在(即已tracked), 则跳过.
              {
                continue;
              }

              //if (marked_views.find(view_new) != marked_views.end())
              //{
              //  continue;
              //}
              marked_views.insert(view_new);
              track_set.insert(view_new);	//保存未处理的view
              view_queue.push(view_new);	//加入队列以进行下一次处理.
              image_touched.insert(image_id);	//至此, 已遍历一个邻近图像, 并从中取出与当前View中的同名点一致的一个key, 组成一对邻近图像的View.
            }// for (; itr_neighbor != itr_neighbor_end; ++itr_neighbor) //将一个view的所有相邻匹配都添加至track_set
          }// if (itr_neighbors != image_neighbor_indexer.end())
        }// while (!view_queue.empty()) //
        if (track_set.size() > 1 && image_touched.size() == track_set.size())	//已迭代遍历当前view的所有相关的view,
        {
          Track track;
          auto itr_track = track_set.begin();
          auto itr_track_end = track_set.end();
          bool is_error = false;
          for (; itr_track != itr_track_end; ++itr_track)	//遍历track_set, 
          {
            if (error_views.find(*itr_track) != error_views.end())	//若在 error_views 中存在则跳过
            {
              is_error = true;
              break;
            }
            track.push_back(*itr_track);	//保存至track(track的数据结构与view相同)
          }

          if (!is_error)
          {
            tracks_loose.push_back(track);	//有效的track保存至tracks_loose
          }
          else
          {
            itr_track = track_set.begin();
            itr_track_end = track_set.end();
            for (; itr_track != itr_track_end; ++itr_track)
            {
              error_views.insert(*itr_track);
            }
          }
        }
        else
        {
          auto itr_track = track_set.begin();
          auto itr_track_end = track_set.end();
          for (; itr_track != itr_track_end; ++itr_track)	//track_set数与image_touched数不符的认为是无效的view或track
          {
            error_views.insert(*itr_track);
          }
        }
      }// for (; itr_key_pair != itr_key_pair_end; ++itr_key_pair)
    } //for (; itr_image_pair_normalized != itr_image_pair_normalized_end;
      //     ++itr_image_pair_normalized)

    for (size_t i = 0; i < tracks_loose.size(); i++)	//检查tracks_loose中的每个成员track中的view是否存在于error_view
    {
      bool has_error = 0;
      for (size_t j = 0; j < tracks_loose[i].size(); j++)
      {
        if (error_views.find(tracks_loose[i][j]) != error_views.end())
        {
          has_error = true;
          break;
        }
      }
      if (has_error)
      {
        tracks_loose[i].clear();
      }
    }

    for (size_t i = 0; i < tracks_loose.size(); i++)
    {
      if (tracks_loose[i].size() > 1)
      {
        tracks.push_back(tracks_loose[i]);
      }
    }

    return 0;
#endif
  }

  Err operator() (const TrackContainer& tracks,
                  MatchContainer& matches) const	//从每个track的同名点列表中构造相邻的"匹配对"
  {
    matches.clear();
    size_t number_of_tracks = tracks.size();
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = tracks[i].size();
      for (size_t j = 0; j < number_of_views - 1; j++)
      {
        size_t image0_id = tracks[i][j].first;
        size_t key0_id = tracks[i][j].second;
        size_t image1_id = tracks[i][j + 1].first;
        size_t key1_id = tracks[i][j + 1].second;
        matches[std::make_pair(image0_id, image1_id)].push_back(
          std::make_pair(key0_id, key1_id));
      }
    }
    return 0;
  }
};

}
}

#endif
