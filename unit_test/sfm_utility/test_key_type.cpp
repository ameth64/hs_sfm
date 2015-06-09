#include <fstream>

#include <gtest/gtest.h>

#include <cereal/archives/portable_binary.hpp>

#include "hs_test_utility/test_toolkit/test_serialization.hpp"

#include "hs_sfm/sfm_utility/key_type.hpp"

namespace
{

TEST(TestKeyType, SerializationTest)
{
  typedef double Scalar;
  typedef hs::sfm::ImageKeys<Scalar> Keyset;
  typedef Keyset::Key Key;

  Keyset keyset;
  keyset.AddKey(Key(2.0, 3.0));
  keyset.AddKey(Key(4.0, 5.0));
  keyset.AddKey(Key(6.0, 7.0));

  hs::test::TestLocalSerialization<Keyset> keyset_tester;
  ASSERT_EQ(true, keyset_tester(keyset,
                                "keyset.bin"));
}

}
