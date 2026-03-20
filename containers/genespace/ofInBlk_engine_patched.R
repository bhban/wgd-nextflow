ofInBlk_engine <- function (gsParam, genome1, genome2, overwrite) 
{
    query <- target <- ofID1 <- ofID2 <- blkID <- rid <- uid1 <- uid2 <- regID <- ofID1 <- ofID2 <- sameHog <- hog1 <- hog2 <- sameInblkOG <- hasSelf <- isSyntenic <- noAnchor <- isArrayRep1 <- isArrayRep2 <- NULL
    blNames <- c("ofID1", "ofID2", "pid", "length", "mismatches", 
        "gapopenings", "queryStart", "queryEnd", "subjectStart", 
        "subjectEnd", "Evalue", "bitScore", "rid")
    blNamesR <- c("ofID2", "ofID1", "pid", "length", "mismatches", 
        "gapopenings", "subjectStart", "subjectEnd", "queryStart", 
        "queryEnd", "Evalue", "bitScore", "rid")
    nCores <- gsParam$params$nCores
    md <- data.table(gsParam$synteny$blast)
    x <- subset(md, query == genome1 & target == genome2)
    if (nrow(x) == 0) {
        tmp <- genome2
        genome2 <- genome1
        genome1 <- tmp
    }
    if (nrow(x) > 1) 
        stop(sprintf("genome1: %s and genome2: %s are not unique in gsParam\n", 
            genome1, genome2))
    allBlast <- subset(read_allBlast(x$allBlast), !is.na(regID) & 
        !noAnchor & isArrayRep1 & isArrayRep2)
    sb01 <- data.table(allBlast)
    dontRerun <- any(!is.na(allBlast$sameInblkOG))
    if (dontRerun) 
        dontRerun <- all(!is.na(allBlast$sameInblkOG)) & !overwrite
    if (dontRerun) {
        cat("found existing run, not rerunning\n")
        return(NULL)
    }
    sb01[, `:=`(hasSelf, any(ofID1 %in% ofID2)), by = "regID"]
    onlySelf <- sum(!sb01$hasSelf) < gsParam$params$blkSize
    sb01 <- subset(sb01, !hasSelf)
    sb01[, `:=`(ofID1 = sprintf("%s_g1", ofID1), ofID2 = sprintf("%s_g2", 
        ofID2))]
    if (onlySelf) {
        cat("no non-self blocks\n")
        return(NULL)
    }
    sb01[, `:=`(rid, sprintf("reg%s", as.numeric(as.factor(regID))))]
    targetGenome <- sb01$genome2[1]
    queryGenome <- sb01$genome1[1]
    sb01[, `:=`(uid1 = uniqueN(ofID1), uid2 = uniqueN(ofID2)), 
        by = "rid"]
    sb01md <- sb01[, list(n = (uid1[1] + uid2[1])/2), by = "rid"]
    propPass <- sum(sb01md$n[sb01md$n >= 40])/sum(sb01md$n)
    if (propPass < 0.5) {
        cat(sprintf("<50%% (%s%%) of syn. hits in regions with >=40 genes, not running\n", 
            round(propPass, 3) * 100))
        return(NULL)
    }
    sb01 <- subset(sb01, uid1 >= 40 & uid2 >= 40)
    di1 <- sb01$ofID1
    names(di1) <- sb01$id1
    di2 <- sb01$ofID2
    names(di2) <- sb01$id2
    peps0 <- read_aaFasta(file.path(gsParam$paths$peptide, sprintf("%s.fa", 
        queryGenome)))
    peps1 <- read_aaFasta(file.path(gsParam$paths$peptide, sprintf("%s.fa", 
        targetGenome)))
    peps0 <- peps0[unique(sb01$id1)]
    peps1 <- peps1[unique(sb01$id2)]
    names(peps0) <- di1[names(peps0)]
    names(peps1) <- di2[names(peps1)]
    bl01 <- sb01[, blNames, with = F]
    bl10 <- sb01[, blNamesR, with = F]
    setnames(bl10, blNames)
    #if (targetGenome == queryGenome) {
    #    bl00 <- data.table(bl01)
    #    bl11 <- data.table(bl10)
    #}
    #else {
        p0md <- subset(md, query == target & query == queryGenome)
        p1md <- subset(md, query == target & query == targetGenome)
        bl00 <- fread(p0md$allBlast, showProgress = F, na.strings = c("NA", 
            ""), select = blNames[1:12])
        bl00[, `:=`(ofID1 = sprintf("%s_g1", ofID1), ofID2 = sprintf("%s_g1", 
            ofID2))]
        bl11 <- fread(p1md$allBlast, showProgress = F, na.strings = c("NA", 
            ""), select = blNames[1:12])
        bl11[, `:=`(ofID1 = sprintf("%s_g2", ofID1), ofID2 = sprintf("%s_g2", 
            ofID2))]
    #}
    tmpDir <- gsParam$paths$tmp
    if (dir.exists(tmpDir)) 
        unlink(tmpDir, recursive = T)
    dir.create(tmpDir)
    on.exit(expr = unlink(list.files(tmpDir, full.names = T), 
        recursive = T))
    blkIDs <- unique(bl01$rid)
    tmpDirs <- sapply(blkIDs, function(j) {
        pth <- file.path(tmpDir, j)
        if (dir.exists(pth)) 
            unlink(pth, recursive = T)
        dir.create(pth)
        dir.create(file.path(pth, "pep"))
        return(pth)
    })
    spl01 <- split(sb01, by = "rid")
    for (j in blkIDs) {
        y <- spl01[[j]]
        p0 <- peps0[unique(y$ofID1)]
        p1 <- peps1[unique(y$ofID2)]
        writeXStringSet(p0, filepath = file.path(tmpDirs[j], 
            "pep", "g0.fa"))
        writeXStringSet(p1, filepath = file.path(tmpDirs[j], 
            "pep", "g1.fa"))
    }
    ofcall <- gsParam$shellCalls$orthofinder
    ofDirs <- rbindlist(mclapply(names(spl01), mc.cores = nCores, 
        function(j) {
            pdir <- file.path(tmpDirs[j], "pep")
            odir <- file.path(tmpDirs[j], "orthofinder")
            if (dir.exists(odir)) 
                unlink(odir, recursive = T)
            ofComm <- sprintf("-f %s -op -o %s", pdir, odir)
            outp <- system2(ofcall, ofComm, stdout = TRUE, stderr = TRUE)
            fp <- list.files(path = odir, pattern = "SequenceIDs.txt", 
                full.names = T, recursive = T)
            out <- data.table(regID = j, path = dirname(fp))
            return(out)
        }))
    write_blast <- function(blastHits, filepath) {
        fwrite(blastHits, file = filepath, sep = "\t", quote = F, 
            row.names = F, col.names = F, showProgress = F)
    }
    read_hogog <- function(path, genomeIDs, allowInBlkOGs) {
        HOG <- list.files(path, pattern = "N0.tsv", recursive = T, 
            full.names = T)
        if (length(HOG) == 1) {
            ogOut <- parse_hogs(filepath = HOG)
            setnames(ogOut, 1, "ogID")
            return(ogOut)
        }
        else {
            if (allowInBlkOGs) {
                OG <- list.files(path, pattern = "Orthogroups.tsv", 
                  recursive = T, full.names = T)
                if (length(OG) == 1) {
                  ogOut <- parse_ogs(filepath = OG, genomeIDs = genomeIDs)
                }
                else {
                  ogOut <- NULL
                }
            }
        }
        return(ogOut)
    }
    inblkHOGs <- rbindlist(mclapply(1:nrow(ofDirs), mc.cores = nCores, 
        function(k) {
            j <- ofDirs$regID[k]
            regj <- ofDirs$path[k]
            sidf <- file.path(regj, "SequenceIDs.txt")
            if (file.exists(sidf)) {
                sids <- read_orthofinderSequenceIDs(sidf)
                idv <- sids$ofID
                names(idv) <- sids$id
                y <- data.table(spl01[[j]])
                u1 <- unique(y$ofID1)
                u2 <- unique(y$ofID2)
                y01 <- subset(bl01, ofID1 %in% u1 & ofID2 %in% 
                  u2)
                y10 <- subset(bl10, ofID1 %in% u2 & ofID2 %in% 
                  u1)
                y00 <- subset(bl00, ofID1 %in% u1 & ofID2 %in% 
                  u1)
                y11 <- subset(bl11, ofID1 %in% u2 & ofID2 %in% 
                  u2)
                y01[, `:=`(ofID1 = idv[ofID1], ofID2 = idv[ofID2])]
                y10[, `:=`(ofID1 = idv[ofID1], ofID2 = idv[ofID2])]
                y00[, `:=`(ofID1 = idv[ofID1], ofID2 = idv[ofID2])]
                y11[, `:=`(ofID1 = idv[ofID1], ofID2 = idv[ofID2])]
                write_blast(y01, file.path(regj, "Blast0_1.txt.gz"))
                write_blast(y10, file.path(regj, "Blast1_0.txt.gz"))
                write_blast(y00, file.path(regj, "Blast0_0.txt.gz"))
                write_blast(y11, file.path(regj, "Blast1_1.txt.gz"))
                comm <- sprintf("-b %s -t 1 -a 1 -X", regj)
                outp <- system2(ofcall, comm, stdout = TRUE, 
                  stderr = TRUE)
                ogout <- read_hogog(path = regj, genomeIDs = c("g0", 
                  "g1"), allowInBlkOGs = TRUE)
                if (!is.null(ogout)) 
                  ogout[, `:=`(regID, j)]
            }
            else {
                ogout <- NULL
            }
            return(ogout)
        }))
    splhog <- split(inblkHOGs, by = "regID")
    hogout <- rbindlist(lapply(names(splhog), function(j) {
        y <- spl01[[j]]
        z <- splhog[[j]]
        if (length(z) == 4) {
            hogv <- z$ogID
            names(hogv) <- z$id
            y[, `:=`(hog1 = hogv[ofID1], hog2 = hogv[ofID2])]
            y[, `:=`(sameHog, hog1 == hog2)]
            return(subset(y, sameHog)[, c("ofID1", "ofID2")])
        }
        else {
            return(NULL)
        }
    }))
    u <- with(hogout, paste(gsub("_g1$", "", ofID1), gsub("_g2$", 
        "", ofID2)))
    allBlast[, `:=`(sameInblkOG, paste(ofID1, ofID2) %in% u)]
    out <- subset(allBlast, sameInblkOG)
    return(out)
}
