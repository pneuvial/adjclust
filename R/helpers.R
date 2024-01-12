matL <- function(mat, h) {
  p <- ncol(mat)
  x <- rep(1:p, each = h+1)
  y <- x + rep(0:h, p)
  coord <- cbind(x, y)
  out <- rep(0, p * (h+1))
  out[y <= p] <- mat[coord[y <= p, ]]
  out <- matrix(out, ncol = h+1, byrow = TRUE)
  out[ ,-1] <- 2*out[ ,-1]
  return(out)
}

matR <- function(mat, h) {
  if ("dgCMatrix" %in% class(mat))
    mat <- t(mat)
  p <- ncol(mat)
  x <- rep(p:1, each = h+1)
  y <- x - rep(0:h, p)
  coord <- cbind(x, y)
  out <- rep(0, p * (h+1))
  out[y >= 1] <- mat[coord[y >= 1, ]]
  out <- matrix(out, ncol = h+1, byrow = TRUE)
  out[ ,-1] <- 2*out[ ,-1]
  return(out)
}

.memberDend <- getFromNamespace(".memberDend", "stats")
.validity.hclust <- getFromNamespace(".validity.hclust", "stats")
plotNodeLimit <- getFromNamespace("plotNodeLimit", "stats")
.midDend <- getFromNamespace(".midDend", "stats")

checkCondition <- function(mat) {# used only with symmetric matrices to check condition
  tmp <- sweep(-2*mat, 1, diag(mat), "+")
  tmp <- sweep(tmp, 2, diag(mat), "+")
  tmp <- tmp[upper.tri(tmp)]
  if (any(tmp < 0)) {
    return(-min(tmp) * 1.01)
  } else return(NULL)
}

#' @importFrom utils getFromNamespace
alt.as.dendrogram <- function(object, hang = -1, check = TRUE, ...) {# used for adding the branching order as a node attribute while turning "hclust" objects into "dendrograms"
  
  t <- 1
  nolabels <- is.null(object$labels)
  merge <- object$merge
  if (check && !isTRUE(msg <- .validity.hclust(object, merge, 
                                               order = nolabels))) 
    stop(msg)
  if (nolabels) 
    object$labels <- seq_along(object$order)
  z <- list()
  nMerge <- length(oHgt <- object$height)
  hMax <- oHgt[nMerge]
  for (k in 1L:nMerge) {
    x <- merge[k, ]
    if (any(neg <- x < 0)) 
      h0 <- if (hang < 0) 
        0
    else max(0, oHgt[k] - hang * hMax)
    if (all(neg)) {
      zk <- as.list(-x)
      attr(zk, "members") <- 2L
      attr(zk, "midpoint") <- 0.5
      attr(zk, "nodePar") <- list(pch=paste(t,sep=""))
      objlabels <- object$labels[-x]
      attr(zk[[1L]], "label") <- objlabels[1L]
      attr(zk[[2L]], "label") <- objlabels[2L]
      attr(zk[[1L]], "members") <- attr(zk[[2L]], "members") <- 1L
      attr(zk[[1L]], "height") <- attr(zk[[2L]], "height") <- h0
      attr(zk[[1L]], "leaf") <- attr(zk[[2L]], "leaf") <- TRUE
      attr(zk[[1L]], "nodePar") <- attr(zk[[2L]], "nodePar") <- list(pch=NA)
      t <- t+1
    }
    else if (any(neg)) {
      X <- as.character(x)
      isL <- x[1L] < 0
      if (isL) {zk <- list(-x[1L], z[[X[2L]]])
      attr(zk[[1L]], "nodePar") <- list(pch=NA)}
      else {zk <- list(z[[X[1L]]], -x[2L])
      attr(zk[[2L]], "nodePar") <- list(pch=NA)}
      attr(zk, "members") <- attr(z[[X[1 + isL]]], "members") + 
        1L
      attr(zk, "midpoint") <- (.memberDend(zk[[1L]]) + 
                                 attr(z[[X[1 + isL]]], "midpoint"))/2
      attr(zk, "nodePar") <- list(pch=paste(t,sep=""))
      attr(zk[[2 - isL]], "members") <- 1L
      attr(zk[[2 - isL]], "height") <- h0
      attr(zk[[2 - isL]], "label") <- object$labels[-x[2 - 
                                                         isL]]
      attr(zk[[2 - isL]], "leaf") <- TRUE
      z[[X[1 + isL]]] <- NULL
      t <- t+1
    }
    else {
      x <- as.character(x)
      zk <- list(z[[x[1L]]], z[[x[2L]]])
      attr(zk, "members") <- attr(z[[x[1L]]], "members") + 
        attr(z[[x[2L]]], "members")
      attr(zk, "midpoint") <- (attr(z[[x[1L]]], "members") + 
                                 attr(z[[x[1L]]], "midpoint") + attr(z[[x[2L]]], 
                                                                     "midpoint"))/2
      attr(zk, "nodePar") <- list(pch=paste(t,sep=""))
      z[[x[1L]]] <- z[[x[2L]]] <- NULL
      t <- t+1
    }
    attr(zk, "height") <- oHgt[k]
    z[[as.character(k)]] <- zk
  }
  structure(z[[as.character(k)]], class = "dendrogram")
}

#' @importFrom stats is.leaf
#' @importFrom graphics strheight strwidth text polygon rect
#' @importFrom utils str
#' @importFrom graphics par
alt.plotNode <- function(x1, x2, subtree, type, center, leaflab, dLeaf, nodePar,
                         edgePar, horiz = FALSE) {#modified plotNode function to be able to print a string instead of a unique symbol at each node
  
  inner <- !is.leaf(subtree) && x1 != x2
  yTop <- attr(subtree, "height")
  bx <- plotNodeLimit(x1, x2, subtree, center)
  xTop <- bx$x
  hasP <- !is.null(nPar <- attr(subtree, "nodePar"))
  if (!hasP) 
    nPar <- nodePar
  if (getOption("verbose")) {
    cat(if (inner) 
      "inner node"
      else "leaf", ":")
    if (!is.null(nPar)) {
      cat(" with node pars\n")
      str(nPar)
    }
    cat(if (inner) 
      paste(" height", formatC(yTop), "; "), "(x1,x2)= (", 
      formatC(x1, width = 4), ",", formatC(x2, width = 4), 
      ")", "--> xTop=", formatC(xTop, width = 8), "\n", 
      sep = "")
  }
  Xtract <- function(nam, L, default, indx) rep(if (nam %in% 
                                                    names(L)) L[[nam]] else default, length.out = indx)[indx]
  asTxt <- function(x) if (is.character(x) || is.expression(x) || 
                           is.null(x)) 
    x
  else as.character(x)
  i <- if (inner || hasP) 
    1
  else 2
  if (!is.null(nPar)) {
    pch <- Xtract("pch", nPar, default = 1L:2, i)
    cex <- Xtract("cex", nPar, default = c(1, 1), i)
    bg <- Xtract("bg", nPar, default = par("bg"), i)
    text(if (horiz)
      cbind(yTop, xTop)
      else cbind(xTop, yTop), labels = pch, adj = c(0.5,-0.75), bg = bg, col = 'blue',
      cex = 0.75*cex, xpd=NA)
  }
  if (leaflab == "textlike") 
    p.col <- Xtract("p.col", nPar, default = "white", i)
  lab.col <- Xtract("lab.col", nPar, default = par("col"), 
                    i)
  lab.cex <- Xtract("lab.cex", nPar, default = c(1, 1), i)
  lab.font <- Xtract("lab.font", nPar, default = par("font"), 
                     i)
  lab.xpd <- Xtract("xpd", nPar, default = c(TRUE, TRUE), i)
  if (is.leaf(subtree)) {
    if (leaflab == "perpendicular") {
      if (horiz) {
        X <- yTop + dLeaf * lab.cex
        Y <- xTop
        srt <- 0
        adj <- c(0, 0.5)
      }
      else {
        Y <- yTop - dLeaf * lab.cex
        X <- xTop
        srt <- 90
        adj <- 1
      }
      nodeText <- asTxt(attr(subtree, "label"))
      text(X, Y, nodeText, xpd = lab.xpd, srt = srt, adj = adj, 
           cex = lab.cex, col = lab.col, font = lab.font)
    }
  }
  else if (inner) {
    segmentsHV <- function(x0, y0, x1, y1) {
      if (horiz) 
        segments(y0, x0, y1, x1, col = col, lty = lty, 
                 lwd = lwd)
      else segments(x0, y0, x1, y1, col = col, lty = lty, 
                    lwd = lwd)
    }
    for (k in seq_along(subtree)) {
      child <- subtree[[k]]
      yBot <- attr(child, "height")
      if (getOption("verbose")) 
        cat("ch.", k, "@ h=", yBot, "; ")
      if (is.null(yBot)) 
        yBot <- 0
      xBot <- if (center) 
        mean(bx$limit[k:(k + 1)])
      else bx$limit[k] + .midDend(child)
      hasE <- !is.null(ePar <- attr(child, "edgePar"))
      if (!hasE) 
        ePar <- edgePar
      i <- if (!is.leaf(child) || hasE) 
        1
      else 2
      col <- Xtract("col", ePar, default = par("col"), 
                    i)
      lty <- Xtract("lty", ePar, default = par("lty"), 
                    i)
      lwd <- Xtract("lwd", ePar, default = par("lwd"), 
                    i)
      if (type == "triangle") {
        segmentsHV(xTop, yTop, xBot, yBot)
      }
      else {
        segmentsHV(xTop, yTop, xBot, yTop)
        segmentsHV(xBot, yTop, xBot, yBot)
      }
      vln <- NULL
      if (is.leaf(child) && leaflab == "textlike") {
        nodeText <- asTxt(attr(child, "label"))
        # cat("nodeText 2 vaut : ")
        # print(nodeText)
        if (getOption("verbose")) 
          cat("-- with \"label\"", format(nodeText))
        hln <- 0.6 * strwidth(nodeText, cex = lab.cex)/2
        vln <- 1.5 * strheight(nodeText, cex = lab.cex)/2
        rect(xBot - hln, yBot, xBot + hln, yBot + 2 * 
               vln, col = p.col)
        text(xBot, yBot + vln, nodeText, xpd = lab.xpd, srt=120, 
             cex = lab.cex, col = lab.col, font = lab.font)
      }
      if (!is.null(attr(child, "edgetext"))) {
        edgeText <- asTxt(attr(child, "edgetext"))
        if (getOption("verbose")) 
          cat("-- with \"edgetext\"", format(edgeText))
        if (!is.null(vln)) {
          mx <- if (type == "triangle") 
            (xTop + xBot + ((xTop - xBot)/(yTop - yBot)) * 
               vln)/2
          else xBot
          my <- (yTop + yBot + 2 * vln)/2
        }
        else {
          mx <- if (type == "triangle") 
            (xTop + xBot)/2
          else xBot
          my <- (yTop + yBot)/2
        }
        p.col <- Xtract("p.col", ePar, default = "white", 
                        i)
        p.border <- Xtract("p.border", ePar, default = par("fg"), 
                           i)
        p.lwd <- Xtract("p.lwd", ePar, default = lwd, 
                        i)
        p.lty <- Xtract("p.lty", ePar, default = lty, 
                        i)
        t.col <- Xtract("t.col", ePar, default = col, 
                        i)
        t.cex <- Xtract("t.cex", ePar, default = 1, i)
        t.font <- Xtract("t.font", ePar, default = par("font"), 
                         i)
        vlm <- strheight(c(edgeText, "h"), cex = t.cex)/2
        hlm <- strwidth(c(edgeText, "m"), cex = t.cex)/2
        hl3 <- c(hlm[1L], hlm[1L] + hlm[2L], hlm[1L])
        if (horiz) {
          polygon(my + c(-hl3, hl3), mx + sum(vlm) * 
                    c(-1L:1L, 1L:-1L), col = p.col, border = p.border, 
                  lty = p.lty, lwd = p.lwd)
          text(my, mx, edgeText, cex = t.cex, col = t.col, 
               font = t.font)
        }
        else {
          polygon(mx + c(-hl3, hl3), my + sum(vlm) * 
                    c(-1L:1L, 1L:-1L), col = p.col, border = p.border, 
                  lty = p.lty, lwd = p.lwd)
          text(mx, my, edgeText, cex = t.cex, col = t.col, 
               font = t.font)
        }
      }
      alt.plotNode(bx$limit[k], bx$limit[k + 1], subtree = child, 
                   type, center, leaflab, dLeaf, nodePar, edgePar, 
                   horiz)
    }
  }
  invisible()
}

#' @importFrom stats is.leaf
#' @importFrom grDevices dev.flush dev.hold
#' @importFrom graphics strheight strwidth text
alt.plot <- function(x, type = c("rectangle", "triangle"), center = FALSE, 
                     edge.root = is.leaf(x) || !is.null(attr(x, "edgetext")), 
                     nodePar = NULL, edgePar = list(), leaflab = c("perpendicular", 
                                                                   "textlike", "none"), dLeaf = NULL, xlab = "", ylab = "", 
                     xaxt = "n", yaxt = "s", horiz = FALSE, frame.plot = FALSE, 
                     xlim, ylim, ...) {#To use alt.plotNode instead of plotNode
  
  type <- match.arg(type)
  leaflab <- match.arg(leaflab)
  hgt <- attr(x, "height")
  if (edge.root && is.logical(edge.root)) 
    edge.root <- 0.0625 * if (is.leaf(x)) 
      1
  else hgt
  mem.x <- .memberDend(x)
  yTop <- hgt + edge.root
  if (center) {
    x1 <- 0.5
    x2 <- mem.x + 0.5
  }
  else {
    x1 <- 1
    x2 <- mem.x
  }
  xl. <- c(x1 - 1/2, x2 + 1/2)
  yl. <- c(0, yTop)
  if (horiz) {
    tmp <- xl.
    xl. <- rev(yl.)
    yl. <- tmp
    tmp <- xaxt
    xaxt <- yaxt
    yaxt <- tmp
  }
  if (missing(xlim) || is.null(xlim)) 
    xlim <- xl.
  if (missing(ylim) || is.null(ylim)) 
    ylim <- yl.
  dev.hold()
  on.exit(dev.flush())
  plot(0, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, 
       ylab = ylab, xaxt = xaxt, yaxt = yaxt, frame.plot = frame.plot, 
       ...)
  if (is.null(dLeaf)) 
    dLeaf <- 0.75 * (if (horiz) 
      strwidth("w")
      else strheight("x"))
  if (edge.root) {
    x0 <- plotNodeLimit(x1, x2, x, center)$x
    if (horiz) 
      segments(hgt, x0, yTop, x0)
    else segments(x0, hgt, x0, yTop)
    if (!is.null(et <- attr(x, "edgetext"))) {
      my <- mean(hgt, yTop)
      if (horiz) 
        text(my, x0, et)
      else text(x0, my, et)
    }
  }
  alt.plotNode(x1, x2, x, type = type, center = center, leaflab = leaflab,
               dLeaf = dLeaf, nodePar = nodePar, edgePar = edgePar, 
               horiz = horiz)
}

# update call to replace e.g. 'run.adjclust' by 'adjclust'
update_call <- function(x, name_to) {
  lst <- as.list(x)
  lst[[1]] <- as.symbol(name_to)
  as.call(lst)
}

# compute dendrogram heights depending on the mode (used in plotSim and chac)
dendro_for_mode <- function(x, 
                            mode = c("standard", "corrected", "total-disp", 
                                     "within-disp", "average-disp")) {
  
  if ((x$method == "adjClust-corrected") & (mode != "standard")) {
    stop("Already corrected 'chac' object. 'mode' must be set to 'standard'")
  }
  
  if (mode == "standard") {
    if (any(diff(x$height) < 0)) 
      warning(paste0("\nDetected reversals in dendrogram: ",
                     "mode = 'corrected', 'within-disp' or 'total-disp' might be more relevant."), 
              call. = FALSE)
  } else if (mode == "corrected") {
    x <- correct(x)
  } else if (mode == "total-disp") {
    x$height <- cumsum(x$height)
  } else if (mode == "within-disp" | mode == "average-disp") {
    to_correct <- which((x$merge[ ,1] > 0) | (x$merge[ ,2] > 0))
    for (ind in to_correct) {
      clusters <- x$merge[ind, ]
      clusters <- clusters[clusters > 0]
      x$height[ind] <- x$height[ind] + sum(x$height[clusters])
    }
    if (mode == "average-disp") {
      # search for current cluster size
      out <- sapply((length(x$height) + 1):1, function(num) {
        tmp <- table(table(cutree(as.hclust(x), k = num)))
        res <- rep(0, length(x$height) + 1)
        res[as.numeric(names(tmp))] <- tmp
        return(res)
      })
      out <- apply(out, 1, diff)
      sizes <- unlist(apply(out, 1, function(acol) which(acol > 0)))
      sizes <- as.numeric(sizes)
      x$height <- x$height / sizes
      if (any(diff(x$height) < 0)) 
        warning(paste0("\nDetected reversals in dendrogram: ",
                       "mode = 'corrected', 'within-disp' or 'total-disp' might be more relevant."))
    }
  }
  
  return(x)
}