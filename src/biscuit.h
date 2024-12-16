/* biscuit.h
 *
 * The MIT License (MIT)
 * Copyright (c) 2023-2024 Jacob.Morrison@vai.org
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
**/

#ifndef BISCUIT_VERSION_H
#define BISCUIT_VERSION_H

// Macros to turn version numbers into string literals
#define TO_STR_EX(n) #n
#define TO_STR(n) TO_STR_EX (n)

#define BISCUIT_VERSION_MAJOR 1
#define BISCUIT_VERSION_MINOR 6
#define BISCUIT_VERSION_PATCH 1
#define BISCUIT_VERSION_DEVEL "-dev"

#define PACKAGE_VERSION \
    TO_STR (BISCUIT_VERSION_MAJOR) \
    "." \
    TO_STR (BISCUIT_VERSION_MINOR) \
    "." \
    TO_STR (BISCUIT_VERSION_PATCH) \
    BISCUIT_VERSION_DEVEL

#endif /* BISCUIT_VERSION_H */
