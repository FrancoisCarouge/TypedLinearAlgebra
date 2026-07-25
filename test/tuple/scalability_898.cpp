/* Typed Linear Algebra
Version 0.2.0
https://github.com/FrancoisCarouge/TypedLinearAlgebra

SPDX-License-Identifier: Unlicense

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <https://unlicense.org> */

#include <tuple>

namespace fcarouge::test {
namespace {
//! @test Explore the scalability of the tuple list.
[[maybe_unused]] const auto test{[] {
  using tuple = std::tuple<                             // 0
      int, int, int, int, int, int, int, int, int, int, // 10
      int, int, int, int, int, int, int, int, int, int, // 20
      int, int, int, int, int, int, int, int, int, int, // 30
      int, int, int, int, int, int, int, int, int, int, // 40
      int, int, int, int, int, int, int, int, int, int, // 50
      int, int, int, int, int, int, int, int, int, int, // 60
      int, int, int, int, int, int, int, int, int, int, // 70
      int, int, int, int, int, int, int, int, int, int, // 80
      int, int, int, int, int, int, int, int, int, int, // 90
      int, int, int, int, int, int, int, int, int, int, // 100
      int, int, int, int, int, int, int, int, int, int, // 110
      int, int, int, int, int, int, int, int, int, int, // 120
      int, int, int, int, int, int, int, int, int, int, // 130
      int, int, int, int, int, int, int, int, int, int, // 140
      int, int, int, int, int, int, int, int, int, int, // 150
      int, int, int, int, int, int, int, int, int, int, // 160
      int, int, int, int, int, int, int, int, int, int, // 170
      int, int, int, int, int, int, int, int, int, int, // 180
      int, int, int, int, int, int, int, int, int, int, // 190
      int, int, int, int, int, int, int, int, int, int, // 200
      int, int, int, int, int, int, int, int, int, int, // 210
      int, int, int, int, int, int, int, int, int, int, // 220
      int, int, int, int, int, int, int, int, int, int, // 230
      int, int, int, int, int, int, int, int, int, int, // 240
      int, int, int, int, int, int, int, int, int, int, // 250
      int, int, int, int, int, int, int, int, int, int, // 260
      int, int, int, int, int, int, int, int, int, int, // 270
      int, int, int, int, int, int, int, int, int, int, // 280
      int, int, int, int, int, int, int, int, int, int, // 290
      int, int, int, int, int, int, int, int, int, int, // 300
      int, int, int, int, int, int, int, int, int, int, // 310
      int, int, int, int, int, int, int, int, int, int, // 320
      int, int, int, int, int, int, int, int, int, int, // 330
      int, int, int, int, int, int, int, int, int, int, // 340
      int, int, int, int, int, int, int, int, int, int, // 350
      int, int, int, int, int, int, int, int, int, int, // 360
      int, int, int, int, int, int, int, int, int, int, // 370
      int, int, int, int, int, int, int, int, int, int, // 380
      int, int, int, int, int, int, int, int, int, int, // 390
      int, int, int, int, int, int, int, int, int, int, // 400
      int, int, int, int, int, int, int, int, int, int, // 410
      int, int, int, int, int, int, int, int, int, int, // 420
      int, int, int, int, int, int, int, int, int, int, // 430
      int, int, int, int, int, int, int, int, int, int, // 440
      int, int, int, int, int, int, int, int, int, int, // 450
      int, int, int, int, int, int, int, int, int, int, // 460
      int, int, int, int, int, int, int, int, int, int, // 470
      int, int, int, int, int, int, int, int, int, int, // 480
      int, int, int, int, int, int, int, int, int, int, // 490
      int, int, int, int, int, int, int, int, int, int, // 500
      int, int, int, int, int, int, int, int, int, int, // 510
      int, int, int, int, int, int, int, int, int, int, // 520
      int, int, int, int, int, int, int, int, int, int, // 530
      int, int, int, int, int, int, int, int, int, int, // 540
      int, int, int, int, int, int, int, int, int, int, // 550
      int, int, int, int, int, int, int, int, int, int, // 560
      int, int, int, int, int, int, int, int, int, int, // 570
      int, int, int, int, int, int, int, int, int, int, // 580
      int, int, int, int, int, int, int, int, int, int, // 590
      int, int, int, int, int, int, int, int, int, int, // 600
      int, int, int, int, int, int, int, int, int, int, // 610
      int, int, int, int, int, int, int, int, int, int, // 620
      int, int, int, int, int, int, int, int, int, int, // 630
      int, int, int, int, int, int, int, int, int, int, // 640
      int, int, int, int, int, int, int, int, int, int, // 650
      int, int, int, int, int, int, int, int, int, int, // 660
      int, int, int, int, int, int, int, int, int, int, // 670
      int, int, int, int, int, int, int, int, int, int, // 680
      int, int, int, int, int, int, int, int, int, int, // 690
      int, int, int, int, int, int, int, int, int, int, // 700
      int, int, int, int, int, int, int, int, int, int, // 710
      int, int, int, int, int, int, int, int, int, int, // 720
      int, int, int, int, int, int, int, int, int, int, // 730
      int, int, int, int, int, int, int, int, int, int, // 740
      int, int, int, int, int, int, int, int, int, int, // 750
      int, int, int, int, int, int, int, int, int, int, // 760
      int, int, int, int, int, int, int, int, int, int, // 770
      int, int, int, int, int, int, int, int, int, int, // 780
      int, int, int, int, int, int, int, int, int, int, // 790
      int, int, int, int, int, int, int, int, int, int, // 800
      int, int, int, int, int, int, int, int, int, int, // 810
      int, int, int, int, int, int, int, int, int, int, // 820
      int, int, int, int, int, int, int, int, int, int, // 830
      int, int, int, int, int, int, int, int, int, int, // 840
      int, int, int, int, int, int, int, int, int, int, // 850
      int, int, int, int, int, int, int, int, int, int, // 860
      int, int, int, int, int, int, int, int, int, int, // 870
      int, int, int, int, int, int, int, int, int, int, // 880
      int, int, int, int, int, int, int, int, int, int, // 890
      int, int, int, int, int, int, int, int            // 898
      >;

  [[maybe_unused]] tuple t{};

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
