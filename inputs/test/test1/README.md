# INFO (TEST MATRIX)

- the number of nonzeros is 21 (written in the matrix market file)
- the sum of the entries in the matrix is 33
- header: `%%MatrixMarket matrix coordinate real symmetric`

```
A = [
  4  -1   0  -1   0   0   0   0   0
 -1   4  -1   0  -1   0   0   0   0
  0  -1   4   0   0  -1   0   0   0
 -1   0   0   4  -1   0  -1   0   0
  0  -1   0  -1   4  -1   0  -1   0
  0   0  -1   0  -1   4   0   0  -1
  0   0   0  -1   0   0   4  -1   0
  0   0   0   0  -1   0  -1   4  -1
  0   0   0   0   0  -1   0  -1   4
]
```

## RESULT

```
[dev@12a87a6e8a5d work]$ ./build/matrix_loader_test
Loading Matrix Market file: "/work/inputs/test/test1/test1.mtx"
Double precision view
  dimensions: 9 x 9
  nonzeros : 33
  sample entries (row, col, value):
    (0, 0, 4)
    (1, 0, -1)
    (3, 0, -1)
    (0, 1, -1)
    (1, 1, 4)
    (2, 1, -1)
    (4, 1, -1)
    (1, 2, -1)
    (2, 2, 4)
    (5, 2, -1)
    ... 23 more entries

DD precision view
  dimensions: 9 x 9
  nonzeros : 33
  sample entries (row, col, value):
    (0, 0, 4.0000000000000000000000000000000e0       )
    (1, 0, -1.0000000000000000000000000000000e0      )
    (3, 0, -1.0000000000000000000000000000000e0      )
    (0, 1, -1.0000000000000000000000000000000e0      )
    (1, 1, 4.0000000000000000000000000000000e0       )
    (2, 1, -1.0000000000000000000000000000000e0      )
    (4, 1, -1.0000000000000000000000000000000e0      )
    (1, 2, -1.0000000000000000000000000000000e0      )
    (2, 2, 4.0000000000000000000000000000000e0       )
    (5, 2, -1.0000000000000000000000000000000e0      )
    ... 23 more entries

```
