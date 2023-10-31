import re
from typing import Optional
from warnings import warn

from hls4ml.backends.fpga.fpga_types import NamedType
from hls4ml.model.layers import Layer, register_layer
from hls4ml.model.optimizer import OptimizerPass, register_pass
from hls4ml.model.types import FixedPrecisionType

re_purge_prefix = re.compile(r'(?<!\w)(?:ap_|ac_)', re.IGNORECASE)
re_parse_fixed = re.compile(r'\s*(u?)fixed<([^>]+)>\s*', re.IGNORECASE)


class FixedPointQuantizer(Layer):
    def initialize(self):
        inp = self.get_input_variable()
        shape = inp.shape
        dims = inp.dim_names
        self.add_output_variable(shape, dims)
        self.set_attr('n_in', self.get_input_variable().size())
        self.overrides = self.attributes['overrides']
        self.fusible = self.attributes['fusible']
        self.SAT, self.RND = self.attributes['SAT'], self.attributes['RND']
        self.mask_kbi = self.attributes.get('mask_kbi', None)


def to_hls4ml_fixed(fixed: str):
    matched = re_parse_fixed.match(re_purge_prefix.sub('', fixed))
    assert matched is not None, f'Cannot parse {fixed}'
    signed = matched.group(1) != 'u'
    b, i, *args = matched.group(2).split(',')
    b, i = int(b), int(i)
    args = [arg.upper() for arg in args]
    new_type = FixedPrecisionType(b, i, signed, *args)
    # For some reason, __class__ is overwritten in hls4ml
    return new_type


def userconf_ifdef(key: str, layer_name: str, model):
    hls_config: dict = model.config.config['HLSConfig']
    layer_confs: dict = hls_config.get('LayerName', None)
    if not layer_confs:
        return False
    layer_conf = layer_confs.get(layer_name, None)
    if not layer_conf:
        return False
    # return key in layer_conf # Ideal case. Not for now.
    if key.endswith('_t') and key != 'table_t':
        # table_t cannot be defined in Precision, for some reason.
        # On the other hand, result_t, weight_t, bias_t, accum_t cannot be decleared explicitly outside Precision, for now.
        # However, still assume that they can be defined explicitly outside Precision.
        precision_conf = layer_conf.get('Precision', None)
        if not precision_conf:
            return key in layer_conf
        return key[:-2] in precision_conf or key in layer_conf

    if key == 'parallelization_factor':
        # Irregular config key name.
        return 'ParallelizationFactor' in layer_conf

    return key in layer_conf


class EnforceProxyModelEmbeddedConfig(OptimizerPass):
    def match(self, node: Layer):
        if not isinstance(node, FixedPointQuantizer):
            return False
        if not node.overrides:
            return False
        return True

    def transform(self, model, node: FixedPointQuantizer):
        if 'layers' not in node.overrides:
            return False

        graph_changed = False
        layers = node.overrides['layers']
        for name, conf in layers.items():
            conf: dict[str, str]
            name: str
            if name not in model.graph:
                # Some layer may be removed by other passes. (e.g. Final flatten layer)
                continue
            target_node: Layer = model.graph[name]
            for k, v in conf.items():
                if userconf_ifdef(k, name, model):
                    warn(
                        f'Config key {k} is defined in hls_config for layer {name} by user. Proxy model config is ignored.',
                        stacklevel=1,
                    )
                    continue

                if k.endswith('_t'):
                    v0: Optional[NamedType] = target_node.get_attr(k)
                    if v0 is None:
                        continue
                    precision = to_hls4ml_fixed(v)
                    v0.precision = precision
                    v0.name = f'{name}_{k[:-2]}'
                    # Well, it turned out that there is yet ANOTHER copy saved in config....
                    # Whyyyy?
                    model.config.layer_name_precision[f'{name}_{k[:-2]}'] = v
                elif k in target_node.attributes.attributes:
                    target_node.set_attr(k, v)
                elif k == 'parallelization_factor':
                    target_node.set_attr(k, int(v))

            if linear_node := model.graph.get(f'{name}_linear'):
                # Proxy model does not assume any extra linear layer.
                # Purge them on sight
                model.remove_node(linear_node)
                graph_changed = True

        return graph_changed


def register_proxy_model():
    register_layer('FixedPointQuantizer', FixedPointQuantizer)
    register_pass('enforce_proxy_model_embedded_config', EnforceProxyModelEmbeddedConfig)