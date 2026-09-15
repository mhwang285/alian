FROM ubuntu:22.04

LABEL maintainer.name="Tucker Hwang"
LABEL maintainer.email="tucker_hwang@berkeley.edu"

ENV LANG=C.UTF-8

WORKDIR /opt

COPY packages entrypoint.sh requirements.txt apps/

RUN apt-get update -qq \
 && ln -sf /usr/share/zoneinfo/US/Pacific /etc/localtime \
 && apt-get -y install --no-install-recommends $(cat apps/packages) wget \
 && rm -rf /var/lib/apt/lists/* \
 && cd /opt/apps \
 && wget https://sourceforge.net/projects/lmod/files/lua-5.1.4.9.tar.bz2 \
 && tar xf lua-5.1.4.9.tar.bz2 \
 && rm -f lua-5.1.4.9.tar.bz2 \
 && cd lua-5.1.4.9 \
 && ./configure --prefix=/opt/apps/lua/5.1.4.9 \
 && make && make install \
 && cd /opt/apps/lua; ln -s 5.1.4.9 lua \
 && ln -s /opt/apps/lua/lua/bin/lua /usr/local/bin \
 && ln -s /opt/apps/lua/lua/bin/luac /usr/local/bin \
 && cd /opt/apps \
 && wget https://sourceforge.net/projects/lmod/files/Lmod-8.7.tar.bz2 \
 && tar xfvj Lmod-8.7.tar.bz2 \
 && rm -f Lmod-8.7.tar.bz2 \
 && cd Lmod-8.7 \
 && ./configure --prefix=/opt/apps \
 && make install \
 && ln -s /opt/apps/lmod/lmod/init/profile /etc/profile.d/z00_lmod.sh \
 && cd /opt \
 && . /etc/profile.d/z00_lmod.sh \
 && git clone --depth 1 --branch main --single-branch https://github.com/matplo/yasp.git \
 && ./yasp/yaspenv.sh \
   "python3 -m pip install --no-cache-dir -r apps/requirements.txt \
 && yasp -mi bundle/hepbase --opt rootspec=ubuntu22.04 pythiaver=8317 n_cores=2 \
 && module use yasp/software/modules \
 && module load yasp bundle/hepbase \
 && git clone --depth 1 --branch main --single-branch https://github.com/matplo/heppyy.git \
#  && sed -i 's@/tmp@/opt/yasp/.workdir/bundle/hepbase/pythia8/default/pythia8311@g' /opt/heppyy/heppyy/cmake/Modules/FindPythia8.cmake \
 && ./heppyy/install_with_yasp.sh"
COPY alian/ /opt/alian/alian
COPY yasp_recipe/ /opt/alian/yasp_recipe
COPY install_with_yasp.sh /opt/alian/
RUN . /etc/profile.d/z00_lmod.sh \
 && ./yasp/yaspenv.sh \
   "module load heppyy bundle/hepbase \
 && ./alian/install_with_yasp.sh"
ENTRYPOINT [ "/opt/apps/entrypoint.sh" ]
CMD [ "/usr/bin/bash", "-l" ]
