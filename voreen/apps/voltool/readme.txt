voltool synth 64 synth64
voltool conv32 synth64.dat synth64_32bit
voltool histogram synth64.dat synth64.hist

voltool dao 64 synth64
voltool vqpack synth64
voltool vq synth64
voltool vqunpack synth64

voltool vqmeasure synth64



histogram .hist
envogram .env
packed _packed
dao _dao
codebook .cb
index .index