# Code used to prepare the data for unit test Get_MNN_Pairs.
# The data can be visualized to check the soundness of the approach.

x <- matrix(c(1,3,1,1,2,3), ncol = 2, byrow = TRUE)
x


y <- matrix(c(2,1,10,9,20,18,19,15), ncol = 2, byrow = TRUE)
y

d <- rbind(
      data.frame(x, batch = "A"),
      data.frame(y, batch = "B")
)

ggplot(d, aes(X1, X2, color = batch)) +
  geom_point() +
  scale_x_continuous(breaks = 1:20) +
  scale_y_continuous(breaks = 1:20) +
  theme(panel.grid = element_line(color = "grey", size = .1))


foo <- Canek:::Get_MNN_Pairs(t(x), t(y), k_Neighbors = 1)
foo

x[foo$Pairs[, 2], ]


y[foo$Pairs[, 1], ]


foo <- Canek:::Get_MNN_Pairs(t(x), t(y), k_Neighbors = 2)
foo

x[foo$Pairs[, 2], ]


y[foo$Pairs[, 1], ]
