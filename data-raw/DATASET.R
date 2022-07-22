set.seed(1234)

Yinyang1000 = Yinyang_data(d = 1000)
usethis::use_data(Yinyang1000, overwrite = TRUE)

Ring1000 = ring_data(d = 1000)
usethis::use_data(Ring1000, overwrite = TRUE)

Mickey1000 = Mickey_data(d = 1000)
usethis::use_data(Mickey1000, overwrite = TRUE)

MixMickey1000 = MixMickey_data(d = 1000)
usethis::use_data(MixMickey1000, overwrite = TRUE)

MixStar1000 = MixStar_data(d = 1000)
usethis::use_data(MixStar1000, overwrite = TRUE)

ManifoldMix1000 = MM_data(d = 1000)
usethis::use_data(ManifoldMix1000, overwrite = TRUE)
